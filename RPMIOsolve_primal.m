function [optvalue,M]=RPMIOsolve_primal(prob,order)
% prob  -- the problem data
% order -- the order k of the relaxation
% optvalue   -- the optimal value of the relaxation f^dual_k


m=prob.Pdim;
n=prob.Xnum;
l=prob.Ynum;
q=prob.Gdim;
tol=10^(-3);

switch prob.type
    case 'convex'
        k1=order;
    case 'nonconvex'
        k1=order;
    case 'sos-convex'
        k1=ceil(degree(prob.theta)/2);
    case 'linear'
        k1=1;
end
k2=order;

%define the moment sequence S
Sdim=nchoosek(n+2*k2,n);
S=sdpvar(m,m,Sdim);

degs1 = deglist(n, 0, k2);
sdegs1=size(degs1,1);

%compute the moment matrix M_k2(S)
MM=sdpvar(sdegs1*m);
for i=1:sdegs1
    for j=1:i 
        index=getindex(degs1(i,:)+degs1(j,:));
        MM((i-1)*m+1: i*m, (j-1)*m+1: j*m)=S(:,:,index);
        MM((j-1)*m+1: j*m ,(i-1)*m+1: i*m)=S(:,:,index);
    end
end

G=prob.G;
%get the exponents of terms in G and the associated coefficient matrices
C_G=[];% coeffcient matrices 
e_G=[];% exponents of terms in G
for i=1:q
    for j=1:i
        [c_G,t_G]=coefficients(G(i,j),prob.X);
        for k=1:length(c_G)
            if isempty(e_G)
                e_G=[degree(t_G(k),prob.X)];
                temp=zeros(q,q);
                temp(i,j)=c_G(k);
                temp(j,i)=c_G(k);
                C_G=[C_G temp];
            else
                [a,b]=ismember(degree(t_G(k),prob.X), e_G,'row');
                if a==1
                    C_G(i,(b-1)*q+j)=c_G(k);
                    C_G(j,(b-1)*q+i)=c_G(k);
                else
                    e_G=[e_G; degree(t_G(k),prob.X)];
                    temp=zeros(q,q);
                    temp(i,j)=c_G(k);
                    temp(j,i)=c_G(k);
                    C_G=[C_G temp];
                end
            end
        end
    end
end

d_G=max(sum(e_G,2));
scg=size(e_G,1);
    
degsgg = deglist(n, 0, floor(k2-d_G/2));
sdegsgg=size(degsgg,1);
% compute the localizing moment matrix M_k2(GS) 
MMgg=sdpvar(sdegsgg*m*q);
for i=1:sdegsgg
    for j=1:i
        temp=zeros(m*q, m*q);
        for k=1:scg
            index=getindex(degsgg(i,:)+degsgg(j,:)+e_G(k,:));
            temp=temp+kron(S(:,:,index),C_G(:,(k-1)*q+1:k*q));
        end
        MMgg((i-1)*m*q+1: i*m*q, (j-1)*m*q+1: j*m*q)=temp;
        MMgg((j-1)*m*q+1: j*m*q, (i-1)*m*q+1: i*m*q)=temp;
    end
end
F=[MM>=0; MMgg>=0];

f=prob.f;
P=prob.P;
d_x_P=0;
d_y_P=0;
for i=1:m
    for j=1:i
        d_x_P=max(d_x_P, degree(replace(P(i,j),prob.Y,1)));
        d_y_P=max(d_y_P, degree(replace(P(i,j),prob.X,1)));
    end
end

k_x=max(d_G,d_x_P);

%compute L_S(P(y,x))
temp=0;
for i=1:m 
    for j=1:m
        [C,T]=coefficients(P(i,j),prob.X);
        sT=length(T);
        T_exp=[];
        for k=1:sT
            index=getindex(degree(T(k),prob.X));
            temp=temp+C(k)*S(i,j,index);
        end
    end
end

sdpvar rho;

s=[];
c=[];
temp2=0;
if ~isempty(prob.theta)
    temp2=0;
    for i=1:length(prob.theta)
        [s1,c1] = polynomial(prob.Y,2*k1-degree(prob.theta(i)));
        s=[s;s1];
        F=[F, sos(s(i))];
        temp2=temp2+s1*prob.theta(i);
        c=[c;c1];
    end
end


F = [F, sos(f-rho-temp-temp2)];
ss=recover(getvariables(S));


ops = sdpsettings('solver','mosek');
ops = sdpsettings(ops,'verbose',0);
solvesos(F,-rho,ops,[c;ss;rho]);

optvalue=value(rho);
disp(['the primal optimal value f^primal_k at order k=', num2str(order), ' is ', num2str(optvalue)]);

M=value(MM);
flag=0;
for i=ceil(k_x/2):k2
    if rank(M(1:m*nchoosek(n+i, n),1:m*nchoosek(n+i, n)),tol)==...
        rank(M(1:m*nchoosek(n+(i-ceil(d_G/2)), n),1:m*nchoosek(n+(i-ceil(d_G/2)), n)),tol)
        flag=1;
        t=i;
        disp(['the rank condition is satisfied at t=', num2str(i), ' with rank being ', ...
            num2str(rank(M(1:m*nchoosek(n+i, n),1:m*nchoosek(n+i, n)),tol))]);
        if strcmp(prob.type, 'nonconvex')
            disp(['we get LOWER bound of f*: ', num2str(optvalue)]);
        elseif  strcmp(prob.type, 'convex')==false
            disp(['the global optimality is numerically certified ']);
        end
        break;
    end
end



 if flag==1
    disp(['the minimizer S^(k,*) admits a representing measure ']);
    MM=M(1:nchoosek(n+t,t)*m, 1:nchoosek(n+t,t)*m);
    [X, W]=extraction(MM, n, m, t);
 end


