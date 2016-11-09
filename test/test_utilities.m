%% Exhaustive Test
tic;
NtestUtil = 1200;
for i=1:NtestUtil
    n=1+floor(rand*100);
    m=1+floor(rand*80);
    
    x1=rand(n,1);
    x2=rand(m,1);
    Q=rand(m,n);
    A=rand(n,n);
    alpha=1.234;
    beta=7.5;
    b=rand(m,1);
    
    %Q(1:10,1:10)
    [alpha_x1tQx2, alpha_Qx1_pb, alpha_Qtx2_px1, zdot] = ...
        testutil(Q, x1, x2, A, alpha, beta, b);
    
    Jtest = alpha* x1'* Q' * x2 + beta;
    Ftest = alpha*Q*x1+b;
    Ktest = alpha*Q'*x2 + x1;
    Wtest = alpha_Qx1_pb'*x2;
    
    % Assertions:
    assert(norm(alpha_Qtx2_px1-Ktest)<1e-7);
    assert(abs(alpha_x1tQx2-Jtest)<1e-7);
    assert(norm(alpha_Qx1_pb-Ftest)<1e-7);
    assert(norm(zdot-Wtest)<1e-7);
    
end
elapsed_time = toc;
disp(['Test Z-UTILITIES ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);

%% Test ZLOWSLV - TRANSPOSE
clear all
n=5;
tic
L = tril(randn(n,n)); for i=1:n, L(i,i)=i; end; x=randn(n,1);
sol=(L')\x;
z=zlowslv(L,x,1);
error0 = norm((z-sol),'inf');
assert(error0<1e-8, ['zlowslv: inaccurate solution. Inf-Error = ' num2str(error0)]);
elapsed_time = toc;
disp(['Test ZLOWSLV-TR  ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);


%% Test ZLOWSLV
clear all
n=5;
tic
L = tril(randn(n,n)); for i=1:n, L(i,i)=i; end; x=randn(n,1);
sol=L\x;
z=zlowslv(L,x);
error0 = norm((z-sol),'inf');
assert(error0<1e-8, ['zlowslv: inaccurate solution. Inf-Error = ' num2str(error0)]);
elapsed_time = toc;
disp(['Test ZLOWSLV-TR  ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);