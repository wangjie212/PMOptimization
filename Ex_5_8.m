function [prob] = P_5_8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set the RPMIO problem data in the following form
%           min_{y\in Y}  f(y)
%                      s.t. P(y,x)>=0, \forall x\in X\subset R^n                       
%where Y={y\in R^l : theta_1(y)>=0, \ldots, theta_s(y)>=0}
%and X={x\in R^n : G(x)>=0} G(x)\in S[x]^m
%
%Please refer the paper:
%Feng Guo and Jie Wang, A Moment-SOS Hierarchy for Robust 
%Polynomial Matrix Inequality Optimization with SOS-Convexity, 
%arXiv:2304.12628
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
%Number of the parameters x
prob.Xnum=2; 
%Number of the variables y 
prob.Ynum=2; 
%dimension of the matrix P
prob.Pdim=2; 
%dimension of the matrix G
prob.Gdim=4; 
%==========================================================================
%Please use x to denote the variables and y to denote the parameters
%DO NOT change
x=sdpvar(prob.Xnum,1);
y=sdpvar(prob.Ynum,1);
prob.X=x;
prob.Y=y;
%==========================================================================
%the objective f(y)
%prob.f=y(1)*y(2); 
prob.f=(y(1)+1)^2+(y(2)-1)^2;
%prob.f=-y(1)+y(2); %lower bounds
%prob.f=y(1)+y(2) %lower bounds
%prob.f=(1-y(1))^2+(1-y(2))^2; 
%==========================================================================
%Semi-infinite constraint matrix P(y,x)
f=1+(x(1)^2*y(1)^2+x(2)^2*y(2)^2-2*x(1)*x(2)*y(1)*y(2));
g=x(2)*y(1)+x(1)*y(2);
prob.P=[f 2*g; 2*g 1];
%==========================================================================
%index set matrix G(x)
prob.G=[1-x(1)^2 x(2) 0 0 ;x(2) 1 0 0 ; 0 0 x(1)*x(2) 0 ; 0 0 0 x(1)^2+x(2)^2-1]; 
%==========================================================================
%Finitely many constraints
prob.theta=[1-y(1)^2-y(2)^2]; %b=1
%==========================================================================
%Type of the problem. 
%Options:   'convex'; 
%           'sos-convex'; 
%           'linear'; 
%           'nonconvex'
prob.type='nonconvex'; 










