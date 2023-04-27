function [X, W]=extraction(M, n, m, t)
%M -- moment matrix satisfies FEC
%n -- number of x variables
%m -- dimension of P(x,y)
%t -- the order of M at which FEC is satisfied
%X -- atoms [x(1)^T x(2)^T ... x(r)^T]
%W -- weights [W_1,\ldots,W_r]^T

tol=10^(-3);

%compute the Cholesky decomposition
[V,S] = svd(M);
S=diag(S);
ind=find(S>=tol);
V=V(:,ind)*diag(sqrt(S(ind)));

%compute the column echolon form U 
%based on Gloptipoly3's CEFroutine 
U=V;
[U1,U2] = size(U);
i = 1; j = 1; basis = [];
while (i <= U2) & (j <= U1)
 [p,k] = max(abs(U(j,i:U2))); k = k+i-1;
 if (p <= tol)
  U(j,i:U2) = zeros(1,U2-i+1,1);
  j = j + 1;
 else
  basis = [basis j];
  U(j:U1,[i k]) = U(j:U1,[k i]);
  found = 0;
  while ~found
   if abs(U(j,i)) < tol*max(abs(U(:,i)))
    j = j + 1;
    found = (j == U1);
   else
    found = 1;
   end;
  end;
  if j <= U1,
   U(j:U1,i) = U(j:U1,i)/U(j,i);
   for k = [1:i-1 i+1:U2]
    U(j:U1,k) = U(j:U1,k) - U(j,k)*U(j:U1,i);
   end
   i = i + 1;
   j = j + 1;
  end
 end
end


degs1 = deglist(n, 0, t);
WX=[kron(ones(size(degs1,1),1), eye(m)) kron(degs1, ones(m,1))];
xtemp=[zeros(n, m) eye(n)];
[a,b]=ismember(WX(basis,:)+xtemp(1,:), WX, 'row');
N=U(b,:);
for i=2:n
 [a,b]=ismember(WX(basis,:)+xtemp(i,:), WX, 'row');
N=[N; U(b,:)];
end
r=rand(n,1);
r=r/sum(r); 
[A,T] = schur(kron(r',eye(size(basis,2)))*N);
X=[];
for i=1:size(A,2)
 X=[X kron(eye(n),A(:,i)')*N*A(:,i)];
end

%remove the duplicate points x and keep distinct ones 
%X=[x(1)^T x(2)^T ... x(r)^T]
i=1;
while(i<size(X, 2))
 for j= size(X, 2):-1:i+1
  if max(abs(X(:,j)-X(:,i)))<tol
   X(:,j)=[];
  end
end
i=i+1;
end
r= size(X,2);

%compute weight W=[W_1,\ldots,W_r]^T
lambda=[];
for i=1:r
 lambda=[lambda prod(X(:,i)'.^degs1,2)];
end
lambda=kron(lambda, eye(m));


%compute the index set of m*r linear independent rows of lambda
U=lambda;
[U1,U2] = size(U);
i = 1; j = 1; lir = [];
while (i <= U2) & (j <= U1)
 [p,k] = max(abs(U(j,i:U2))); k = k+i-1;
 if (p <= tol)
  U(j,i:U2) = zeros(1,U2-i+1,1);
  j = j + 1;
 else
  lir = [lir j];
  U(j:U1,[i k]) = U(j:U1,[k i]);
  found = 0;
  while ~found
   if abs(U(j,i)) < tol*max(abs(U(:,i)))
    j = j + 1;
    found = (j == U1);
   else
    found = 1;
   end;
  end;
  if j <= U1,
   U(j:U1,i) = U(j:U1,i)/U(j,i);
   for k = [1:i-1 i+1:U2]
    U(j:U1,k) = U(j:U1,k) - U(j,k)*U(j:U1,i);
   end
   i = i + 1;
   j = j + 1;
  end
 end
end
%

if size(lir,2)<m*r
 error('need moment matrix of higher order');
else
 L_R=lambda(lir(1:m*r),:);
end

W=inv(L_R)*M(lir(1:m*r),1:m);


disp(['there are ', num2str(r), ' atoms in the extracted measure: ']);
for i=1:r
        disp(['the ', num2str(i), '-th atom is: ']);
        disp(X(:,i));
        disp(['with the ', num2str(i), '-th weight being']);
        disp(W((i-1)*m+1: i*m, :));
    end
end
