clear all
clear classes
load ./archive/testCase1

max_iter = 250;
e_g = 1e-4; % infeasibility
e_V = 1e-5; % sub-optimality

N=90;

% SOLVE PROBLEM
% x+ = Ax + Bu +f
% xmin <= x <= xmax
% umin <= u <= umax
% cmin <= Fx + Gu <= cmax
x0 = 1.5*rand(nx,1);
% Example
sys = LTISystem(A, B, 0.0*randn(nx,1), F, G, [], g, xmin, xmax, umin, umax);

mpc = MPCGPAD(sys);
mpc=mpc.setTerminalCost(dare(A, B, Q, R, []), []);
mpc.monotone=false;
mpc = mpc.setTerminalConstraints(FN,gN). ...
    setCost(Q, R).setHorizon(N);

[nx, nu, nc] = mpc.sys.getDimensions();
mpc.isInfHorizonCost = 1;
mpc = mpc.factor();
alpha = 1/mpc.calculateLipschitz();

tic;
[Ustar, Xstar, Diag] = mpc.control(1.2*x0, alpha, e_g, e_V, max_iter);
assert(Diag.iter<max_iter, 'Diag.iter >= max_iter - Did not converge');

elapsed_time = toc;
disp(['Test GPAD!       ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);