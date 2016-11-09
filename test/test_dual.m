clear all
load ./archive/testCase1.mat
tic;
%TODO Test with linear terms!!!
for i=1:100
    Q=rand(nx,nx); Q=Q'*Q;
    Q=Q/1000;
    R=rand(nu,nu);R=R'*R;
    Pf=dare(A, B, Q, R, []);
    x0=rand(nx,1);
    np = size(F,1)*2 + 2*nx + 2*nu;
    y=rand(np,N);
    yN=rand(nf,1);
    [K, L, D, M, C, s, d, e, Jdual, Xstar, Ustar] = ...
        solve_dual(A, B, [], Q, R, [], [], [], Pf, [], ...
        F, G, -g, g, xmin, xmax, umin, umax, FN, gN, N, x0, y, yN);
    
    Rbar = R  +  B'*Pf*B;
    Rchol = chol(Rbar,'lower');
    
    [d_gpad, Ustar_gpad,Xstar_gpad] = ...
        dcost(A, B, [], Q, R, [], [], [], Pf, [], FN, gN, F, G, -g, g, xmin, ...
        xmax, umin, umax, D, L, M, K, C, s, d, Rchol, y, yN, x0, N);
    assert(norm(Ustar-Ustar_gpad)<1e-6, 'Ustar : WRONG');
    assert(norm(Xstar-Xstar_gpad)<1e-6, 'Xstar : WRONG');
    assert(norm(Jdual-d_gpad)<1e-6, 'Jdual : WRONG');
end

elapsed_time = toc;
disp(['Test DUALCOST    ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);