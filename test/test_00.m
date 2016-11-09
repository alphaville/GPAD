clear all
load ./archive/testCase1.mat
tic;
clear w
clear wN
for ind=1:40
    Q=rand(nx,nx); Q=Q'*Q;
    Q=Q/1000+1e-3*eye(nx);
    R=rand(nu,nu);R=R'*R;
    Pf=dare(A, B, Q, R, []);
    x0=rand(nx,1);
    np = size(F,1);
    y=rand(np,N);
    yN=rand(nf,1);
    [K, L, D, M, C, s, d, e, Jdual, Xstar, Ustar] = ...
        solve_dual(A, B, [], Q, R, [], [], [], Pf, [], ...
        F, G, [], g, [], [], [], [], FN, gN, N, x0, y, yN);
    
    II=speye((N+1)*nx);
   
    AA = [[zeros(nx,nx*N);kron(speye(N),A)] zeros((N+1)*nx,nx)];
    
    % Equality constraints
    Atilde = AA-II;
    Btilde = [sparse(nx,N*nu);kron(speye(N),B)];
    ftilde = [speye(nx);sparse(N*nx,nx)];
    
    % Inequality constraints
    Ftilde = blkdiag(kron(speye(N),F),FN);
    Gtilde = [kron(speye(N),G);zeros(nf,N*nu)];
    Qtilde = blkdiag(kron(speye(N),Q),Pf);
    Rtilde = kron(speye(N),R);
    Htilde = blkdiag(Qtilde,Rtilde);
    gtilde = [kron(ones(N,1),g);gN];
    
    ytilde = [vec(y);yN];
    [zstar, Jstar, cplexFlag, cplexOutput]= ...
        cplexqp(Htilde, [Ftilde Gtilde]'*ytilde, [], [], [Atilde Btilde], -ftilde*x0);
    Jstar = Jstar - gtilde'*ytilde;
    
    diff_dual_cost = abs(Jdual-Jstar)/abs(Jstar);
    assert(diff_dual_cost<1e-6, ['Dual Cost - Difference = ' num2str(diff_dual_cost)]);
    
    assert(norm(zstar(1:(N+1)*nx)-vec(Xstar))/norm(zstar(1:(N+1)*nx))<1e-6);
    %TODO: Also test with MPCGPAD
    
    e_test=calculate_e(FN, [], L, C, K, [], 1, 0, 1, 1, 1, 1, N, y, yN);
    assert(norm(e_test-e)<1e-6, 'e was not calculated properly in C!');
end
elapsed_time = toc;

disp(['Test DUAL$SOLVE  ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);