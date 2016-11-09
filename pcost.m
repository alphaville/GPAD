%PCOST returns the primal cost given by:
%J = V_f(xN) + sum_k [xk' uk']'[Q S'; S R][xk; uk] + [q' r'][xk; uk]
%where Vf is the terminal cost function given by:
%Vf(x) = x' QN x + qN'x
%for a given sequence of states X and a sequence of inputs U.
%
%Syntax:
%J = PCOST(Q, R, S, q, r, QN, qN, N, X, U, classQ, classR)
%
%Input Arguments:
% Q, R, S : Matrices of the quadratic term of the cost. If Q has the form
% Q = alphaQ*eye(nx) (where nx is the number of states) and alphaQ>=0,
% then we may set classQ=1 and let Q=alphaQ. If Q is a diagonal matrix,
% we may set classQ=2 and let Q be a vector with the diagonal elements 
% thereof. Exactly the same holds for the matrix R. The matrix S is nu-by-nx
% and is considered in general not to be a sparse matrix. If S=[], it will
% not be taken into account in the computations.
%
% N : is the prediction horizon of the MPC problem.
%
% q,r,qN: The vectors q and r, if set to [] (empty), they will be completely 
% discregarded. In most cases it is q=[], r=[], qN=[].
%
% QN is always assumed to be a symmetric matrix. Actually, only the elements
% of its upper triangular part are taken into consideration to minimise
% the iterations. 
%
% X : is a sequence of states represented as a nx-by-K matrix where nx is 
% the dimension of the state vector and K is at least N+1 (The elements of
% X in all columns after the N+1-th - if any - will not be used). Normally,
% X should be nx-by-(N+1) and x0=X(:,1), x1=X(:,2), and so on.
%
% U : is a sequence of input values and is represented as a nu-by-K matrix
% where K is at least equal to N.
%
% classQ, classR: These integer values determine the class of the Matrices
% Q and R respectively. In particular, the following mapping is used:
% 0 : dense matrix/no sparsity
% 1 : Q=alpha*I, this is a quite common case that Q is a multiple of the
%     identity matrix I
% 2 : Q=diag(Qd), i.e., Q is a diagonal matrix
% 3 : Q is a sparse matrix (not implemented yet).
% If classQ or classR are not specified, or they are set to [] (empty),
% the default value 0 is used.
%
%See also:
%make

% Built-in function
% See ./lib_gpad/mex_gpad/mx_pcost.c and ./lib_gpad/gpad.c

% File: pcost.m
% Created on: 15 June, 2013
% Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
% Institute: IMT Lucca, Lucca Italy