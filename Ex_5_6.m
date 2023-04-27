function [prob] = Ex_5_5
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
prob.f=1477307/640000-(341/5000)*y(1)^4*y(2)+(51531/640000)*y(2)^6-(204249/320000)*y(2)^5+(1416373/640000)*y(2)^4+(445373/160000)*y(1)^2-(732223/160000)*y(2)^3-(1091/1250)*y(1)^3+(4112773/640000)*y(2)^2+(6837/40000)*y(1)^4-(641/20000)*y(1)^5+(77/10000)*y(1)^6+(11527/2000)*y(1)*y(2)+(53/2500)*y(1)*y(2)^3-(192337/40000)*y(1)^2*y(2)-(119311/40000)*y(1)*y(2)^2+(26759/80000)*y(1)*y(2)^4+(22717/20000)*y(1)^3*y(2)+(254043/80000)*y(1)^2*y(2)^2-(301/20000)*y(1)^5*y(2)+(2143/40000)*y(1)^4*y(2)^2+(1671/20000)*y(1)^3*y(2)^3+(14901/160000)*y(1)^2*y(2)^4-(1399/20000)*y(1)*y(2)^5-(3281/5000)*y(1)^3*y(2)^2-(33009/40000)*y(1)^2*y(2)^3-(1821097/320000)*y(2)-(272477/80000)*y(1);
%==========================================================================
%Semi-infinite constraint matrix P(y,x)
f=1-(x(1)^2*y(1)^2+x(2)^2*y(2)^2-2*x(1)*x(2)*y(1)*y(2));
g=x(2)*y(1)+x(1)*y(2);
prob.P=[f 2*g; 2*g 1];
%==========================================================================
%index set matrix G(x)
prob.G=[1-x(1)^2 x(2) 0 0 ;x(2) 1 0 0 ; 0 0 x(1)*x(2) 0 ; 0 0 0 x(1)^2+x(2)^2-1]; 
%==========================================================================
%Finitely many constraints
prob.theta=[1-y(1)^2-y(2)^2];
%==========================================================================
%Type of the problem. 
%Options:   'convex'; 
%           'sos-convex'; 
%           'linear'; 
prob.type='convex'; 










