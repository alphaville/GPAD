% Test file for calculate_e
% First run:
% make calculate_e

clear all
load ./archive/testCase1.mat
tic;
for i=1:50
    Q=rand(nx,nx); Q=Q'*Q;
    Q=Q/1000;
    R=rand(nu,nu);R=R'*R;
    Pf=dare(A, B, Q, R, []);
    x0=rand(nx,1);
    np = size(F,1);
    y=rand(np,N);
    yN=rand(nf,1);
    [K, L, ~, ~, C, s, ~, e] = ...
        solve_dual(A, B, [], Q, R, [], [], [], Pf, [], ...
        F, G, [], g, [], [], [], [], FN, gN, N, x0, y, yN);
    e_test = calculate_e(FN, [], L, C, K, s, 1, 0, 1, 1, 1, 1, N, y, yN);
    assert(norm(e-e_test)<1e-7);
    
    
    np = size(F,1)*2;
    y=rand(np,N);
    yN=rand(nf,1);
    [K, L, ~, ~, C, s, ~, e] = ...
        solve_dual(A, B, [], Q, R, [], [], [], Pf, [], ...
        F, G, -g, g, [], [], [], [], FN, gN, N, x0, y, yN);
    e_test = calculate_e(FN, [], L, C, K, s, 0, 0, 1, 1, 1, 1, N, y, yN);
    assert(norm(e-e_test)<1e-7);
    
    np = size(F,1)*2+2*nx+2*nu;
    y=rand(np,N);
    yN=rand(nf,1);
    [K, L, D, M, C, s, d, e, Jdual, Xstar, Ustar] = ...
        solve_dual(A, B, [], Q, R, [], [], [], Pf, [], ...
        F, G, -g, g, xmin, xmax, umin, umax, FN, gN, N, x0, y, yN);
    e_test = calculate_e(FN, [], L, C, K, s, 0, 0, 0, 0, 0, 0, N, y, yN);
    assert(norm(e-e_test)<1e-7);

    np = size(F,1);
    f=rand(nx,1);
    y=rand(np,N);
    yN=rand(nf,1);
    [K, L, ~, ~, C, s, ~, e] = ...
    solve_dual(A, B, f, Q, R, [], [], [], Pf, [], ...
    F, G, [], g, [], [], [], [], FN, gN, N, x0, y, yN);
    e_test = calculate_e(FN, [], L, C, K, s, 1, 0, 1, 1, 1, 1, N, y, yN);
    assert(norm(e-e_test)<1e-7);    
end
elapsed_time = toc;
disp(['Test CALCULATE_E ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);

