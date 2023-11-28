function [xx, rho]=RPMIOsolve_dual(prob,order)
% prob  -- the problem data
% order -- the order k of the relaxation
% xx    -- the approximate minimizer s^(k,*)
% rho   -- the optimal value of the relaxation f^primal_k

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

%define the moment sequence s
degs1 = deglist(l, 0, k1);
sdegs1=size(degs1,1);
sdegs2=nchoosek(l+2*k1,l);
mm=sdpvar(sdegs2,1); 


%compute the moment matrix M_k1(s)
MM=sdpvar(sdegs1);

for i=1:sdegs1
    for j=1:i 
        index=getindex(degs1(i,:)+degs1(j,:));
        MM(i,j)=mm(index);
        MM(j,i)=MM(i,j);
    end
end
F = [mm(1)==1, MM>=0];



%compute the localizing moment matrices at theta
d_theta=0;
for i=1:length(prob.theta)
    d_theta=max(d_theta, degree(prob.theta(i)));
    degs0 = deglist(l, 0, k1-ceil(degree(prob.theta(i))/2));
    sdegs0=size(degs0,1);
    MM0=sdpvar(sdegs0);
    [c_theta,t_theta]=coefficients(prob.theta(i),prob.Y);
    sphi=length(t_theta);
    for j=1:sdegs0
        for k=1:j
            temp=0;
 %           index=getindex(degs0(i,:)+degs0(j,:));
            for w=1:sphi
                index=getindex(degs0(j,:)+degs0(k,:)+degree(t_theta(w),prob.Y));
                temp=temp+c_theta(w)*mm(index);
            end
            MM0(j,k)=temp;
         MM0(k,j)=temp;
        end
    end
    F=[F, MM0>=0];
end


P=prob.P;
d_x_P=0;
d_y_P=0;
for i=1:m
    for j=1:i
        d_x_P=max(d_x_P, degree(replace(P(i,j),prob.Y,1)));
        d_y_P=max(d_y_P, degree(replace(P(i,j),prob.X,1)));
    end
end

G=prob.G;
d_G=0;
for i=1:q
    for j=1:i
        d_G=max(d_G, degree(G(i,j)));        
    end
end


for i=1:m
    for j=1:i
        [C,T]=coefficients(P(i,j),prob.Y);
        sT=length(T);
        temp=0;
        for k=1:sT
            index=getindex(degree(T(k),prob.Y));
            temp=temp+C(k)*mm(index);
        end
        P(i,j)=temp;
        P(j,i)=temp;
    end
end



%define two PSD matrices corresponding to the SOS polynomial matrices in
%quadratic module representation of H_s(P(y,x)) in Q^m_k(G)
v0 = monolist(prob.X,k2);
Z0=sdpvar(nchoosek(n+k2, n)*m);
S0=(kron(v0,eye(m)))'*Z0*(kron(v0,eye(m)));

v1 = monolist(prob.X,floor(k2-d_G/2));
Z1=sdpvar(nchoosek(n+floor(k2-d_G/2), n)*m*q);
S1=(kron(v1,eye(m*q)))'*Z1*(kron(v1,eye(m*q)));


SG=[];
for i=1:m
    SGr=[];
    for j=1:m
        SGr=[SGr sum(sum(S1((i-1)*q+1:i*q,(j-1)*q+1:j*q).*G))];
    end
    SG=[SG; SGr];
end


F = [F, Z0>=0, Z1>=0];
for i=1:m
    for j=1:i
        F=[F, coefficients(P(i,j)-S0(i,j)-SG(i,j),prob.X)==0];
    end
end

%formulate the objective in dual SDP relaxation
obj=prob.f;
[c_f,t_f]=coefficients(obj,prob.Y);
s_f=length(t_f);

k_y=max(max(d_theta, d_y_P),degree(obj));

f=0;
for i=1:s_f
    f=f+c_f(i)*mm(getindex(degree(t_f(i),prob.Y)));
end


ops = sdpsettings('solver','mosek');
ops = sdpsettings(ops,'verbose',0);
optimize(F,f,ops);

%%%%%%%%%%%%%%%%%%%%%%%
%the order of the variables for the minimizer is x(n),...,x(1), so reverse
%%%%%%%%%%%%%%%%%%%%%%
xx=[];
for i=l:-1:1
    xx=[xx value(mm(i+1,1))];
end
xx
rho=value(f)


disp(['the primal dual value f^dual_k at order k=', num2str(order), ' is ', num2str(rho)]);
M=value(MM);


if strcmp(prob.type, 'nonconvex')
    eig(M(1:nchoosek(l+ceil(k_y/2), l),1:nchoosek(l+ceil(k_y/2), l)))
    if rank(M(1:nchoosek(l+ceil(k_y/2), l),1:nchoosek(l+ceil(k_y/2), l)),tol)==1
        disp(['the rank of the moment matrix M_{dy/2} is 1, so we get UPPER bound of f*: ', num2str(rho)]);
    else
        disp(['the rank of the moment matrix M_{dy/2} is NOT 1']);
    end
end
if strcmp(prob.type, 'convex')
    d_Theta=degree(prob.theta);
    flag=0;
    for i=ceil(k_y/2):k1
        if rank(M(1:nchoosek(l+i, l),1:nchoosek(l+i, l)),tol)==...
            rank(M(1:nchoosek(l+(i-ceil(d_Theta/2)), l),1:nchoosek(l+(i-ceil(d_Theta/2)), l)),tol)
            flag=1;
            t=i;
            disp(['the rank condition is satisfied at t=', num2str(i),' with rank being ', ...
                num2str(rank(M(1:nchoosek(l+i, l),1:nchoosek(l+i, l)),tol))]);
            break;
        end
    end
 if flag==1
     disp(['the minimizer s^(k,*) admits a representing measure ']);
     M=M(1:nchoosek(l+t,t), 1:nchoosek(l+t,t));
    [X, W]=extraction(M, l, 1, t);
 end
end







