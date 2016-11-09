%ZLOWSLV solves a triangular system. Given a lower triangular matrix
%L and a vector z, zlowslv solves the linear system Lx=z or the system
%L'x=z by forward and backward substitution respectively. Whether the
%transpose should be used is determined by an input argument. In
%particular, zlowslv has the following structure:
%
%Systax:
% x = zlowslv(L,z,mode);
%
%Input arguments:
% L    :  Is a lower triangular matrix. If the matrix provided is not lower
%         triangular, all other elements will be disregarded.
% z    :  A vector
% mode : If mode is equal to 0, then the system Lx=z will be solved,
%        otherwise if mode is equal to 1 the system L'x=z will be solved
%        instead.
%
%Output arguements:
% x    : Solution to the aforementioned system.
%
%See also:
%test_utilities

% Built-in function
% See ./lib_gpad/mex_gpad/mx_zlowslv.c and ./lib_gpad/gpad.c

% Author: Pantelis Sopasakis
% Created on: 30 September, 2013
% Author: Pantelis Sopasakis <pantelis.sopasakis@imtlucca.it>
% Institute: IMT Lucca, Lucca Italy