++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This is the Matlab scrpt for the numerical experiments in the paper

Feng Guo and Jie Wang, A Moment-SOS Hierarchy for Robust Polynomial Matrix Inequality
Optimization with SOS-Convexity, arXiv: 2304.12628

The primal and dual relaxations in the paper are implemented and
solved for the RPMIO problem 
		min_{y\in Y}  f(y)
                      s.t. P(y,x)>=0, \forall x\in X\subset R^n                       
where Y={y\in R^l : theta_1(y)>=0, \ldots, theta_s(y)>=0}
and X={x\in R^n : G(x)>=0} G(x)\in S[x]^m
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Prerequisites to run the code:

*	Yalmip: to implement the relaxations 
	https://yalmip.github.io/

*	Mosek: to solve the resulting SDP relaxations
	https://www.mosek.com/

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
How to run the code: 

e.g., for the problem in Example 5.6: 
>> prob = Ex_5_6;		%generate the problem data
>> k=3;				%set the order of the relaxation
>> RPMIOsolve_primal(prob,3);	%solve the primal relaxation
>> RPMIOsolve_dual(prob,3);	%solve the dual relaxation
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Problem data files:
Ex_4_14.m:	problem in Example 4.14
Ex_4_18.m:	problem in Example 4.18
Ex_4_21.m:	problem in Example 4.21
Ex_5_6.m:	problem in Example 5.6

Instructions to set the problem data are included in each data file
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Numerical results: Matlab worksheets for each example 
results_Ex_4_14.pdf  
results_Ex_4_18.pdf	
results_Ex_4_21.pdf	
results_Ex_5_6.pdf	
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The script is still at its early stage and rudimentary.
