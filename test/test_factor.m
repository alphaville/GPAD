clear all
clear classes
load ./archive/randomPromlem1
tic;
sys = LTISystem(A, B, [], F, G, [], g, xmin, xmax, umin, umax);
mpc = MPCGPAD(sys);
mpc = mpc.setTerminalConstraints(FN,gN). ...
    setCost(Q, R).setHorizon(N);
mpc = mpc.factor();

P = mpc.QN;
Rbar = R + B'*P*B;
Sbar = B'*P*A;
Rinv = inv(Rbar);
K = -Rinv*Sbar; assert(norm(mpc.K-K)<1e-10);
D = -Rinv*G'; assert(norm(mpc.D-D)<1e-10);
M = -Rinv*B'; assert(norm(mpc.M-M)<1e-10);
L = (A+B*K)'; assert(norm(mpc.L-L)<1e-10);
C = (F+G*K)'; assert(norm(mpc.C-C)<1e-10);
assert(isempty(mpc.s));
assert(isempty(mpc.d));

assert(norm(mpc.RbarChol*mpc.RbarChol'-Rbar)<1e-10);
elapsed_time = toc;
disp(['Test MPC.FACTOR  ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);
