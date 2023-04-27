function [prob] = P_4_20
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
prob.Ynum=1; 
%dimension of the matrix P
prob.Pdim=3; 
%dimension of the matrix G
prob.Gdim=2; 
%==========================================================================
%Please use x to denote the variables and y to denote the parameters
%DO NOT change
x=sdpvar(prob.Xnum,1);
y=sdpvar(prob.Ynum,1);
prob.X=x;
prob.Y=y;
%==========================================================================
%the objective f(y)
prob.f=-y(1); 
%==========================================================================
%Semi-infinite constraint matrix P(y,x)
Q=[1/2^(1/2) -1/3^(1/2) 1/6^(1/2); 0 1/3^(1/2) 2/6^(1/2); 1/2^(1/2) 1/3^(1/2) -1/6^(1/2)];
f1=-x(1)^2-x(2)^2;
f2=-(x(1)+1)^2/4-(x(2)-1)^2/4;
f3=-(x(1)-1)^2/4-(x(2)+1)^2/4;
prob.P=Q*diag([f1,f2,f3])*Q'-y(1)*eye(prob.Pdim);
%==========================================================================
%index set matrix G(x)
prob.G=[1-4*x(1)*x(2) x(1); x(1) 4-x(1)^2-x(2)^2]; 
%==========================================================================
%Finitely many constraints
prob.theta=[];
%==========================================================================
%Type of the problem. 
%Options:   'convex'; 
%           'sos-convex'; 
%           'linear'; 
prob.type='linear'; 










