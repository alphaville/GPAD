%DCOST calculates the dual cost. 
%
%Syntax:
% [d, Ustar, Xstar] = dcost(A, B, f, Q, R, S, q, r, QN, qN, FN, gN, ...
%                           F, G, cmin, cmax, xmin, xmax, umin, umax, ...
%                           D, L, M, K, C, s, d, Rchol, y, yN, x0, N);
%
%Input arguments:
% A, B, f : The matrices that correspond to the affine dynamical system 
%           x+ = Ax + Bu + f
% Q, R, S : Matrices of the quadratic term of the cost. S is allowed to be
%           empty as it is the case quite often
% q,r     : Linear terms of the primal cost; they correspong to the terms
%           q'*x and r'*u in the stage cost. If any of these (or both) are
%           empty, they will be disregarded in the computations (it is
%           equivalent to setting them to zero, nay faster).
% QN      : Matrix of the quadratic part of the terminal cost. In general,
%           the terminal cost is assumed to have the form Vf(x)=x'QN x+qN'x
% qN      : The linear part of the terminal cost. The vector qN is allowed
%           to be empty.
% FN, gN  : The matrix FN and the vector gN define the terminal
%           constraints, i.e., FN xN <= gN
% F,G,cmin, cmax : Users may impose double-sided mixed state-input
%           constraints of the form cmin <= Fx + Gu <= cmax. All these
%           matrices and vectors are allowed to be empty.
% xmin,xmax : Lower and upper boundes on the state vector: xmin <= x <=
%           xmax.
% umin,umax : Lower and upper bounds for the input vector
% D, L, M, K, C, s, d, Rchol : Matrices calculated by the factorisation
%           step offline
% y       : The dual variable for the predicted state sequence up to N-1
% yN      : The dual variable that corresponds to the terminal state
% x0      : The initial (current) state
% N       : The prediction horizon.
%
%
%Output arguments:
%d        : The dual cost
%Xstar    : The predicted sequence of states
%Ustar    : The predicted sequence of inputs
%
%See also:
%pcost, dgrad

% Built-in function
% See ./lib_gpad/mex_gpad/mx_dcost.c and ./lib_gpad/gpad.c

% Author: Pantelis Sopasakis
% Created on: 15 June, 2013
% Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
% Institute: IMT Lucca, Lucca Italy