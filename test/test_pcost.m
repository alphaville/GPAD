%TEST method for PCOST.
tic;
NPcostIter = 500;
for j=1:NPcostIter
    nx=1+floor(120*rand);
    nu=1+floor(80*rand);
    N=10+floor(10*rand);
    X=rand(nx, N+1);
    U=rand(nu, N);
    Q=(5+2*rand)*eye(nx);
    R=(1+rand)*eye(nu);
    if rand>0.5, S=rand(nu, nx); else S=[]; end
    if rand>0.5, q=rand(nx,1); else q=[]; end
    if rand>0.5, r=rand(nu,1); else r=[]; end
    QN=rand(nx,nx);
    QN=QN+QN'; % QN has to be symmetric!
    qN=rand(nx,1);
    J=0.0;
    for i=1:N
        J = J + 0.5*X(:,i)'*Q*X(:,i) + 0.5*U(:,i)'*R*U(:,i);
        if ~isempty(S), J = J + X(:,i)'*S'*U(:,i); end
        if ~isempty(q), J = J + q'*X(:,i); end
        if ~isempty(r), J = J + r'*U(:,i); end
    end
    J = J + 0.5* X(:,N+1)' * QN * X(:,N+1) + qN'*X(:,N+1);
    Jc = pcost(Q, R, S, q, r, QN, qN, N, X, U);
    JcaI = pcost(Q(1), R(1), S, q, r, QN, qN, N, X, U, 1, 1);
    Jcdiag = pcost(diag(Q), diag(R), S, q, r, QN, qN, N, X, U, 2, 2);
    JcdiagQ = pcost(diag(Q), R, S, q, r, QN, qN, N, X, U, 2);
    
    assert(abs(J-Jc)<1e-7, 'PCOST was not calculated accurately! (1)');
    assert(abs(J-JcaI)<1e-7, 'PCOST was not calculated accurately! (2)');
    assert(abs(J-Jcdiag)<1e-7, 'PCOST was not calculated accurately! (3)');
    assert(abs(J-JcdiagQ)<1e-7, 'PCOST was not calculated accurately! (4)');
end
elapsed_time = toc;

disp(['Test PCOST       ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);