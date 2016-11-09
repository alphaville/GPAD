clear all;
tic;
n=100;
for i=1:100
    x=rand(n,1);
    Q=rand(n,n); Q=sqrt(Q*Q'+eye(n));
    L=chol(Q,'lower');
    y=testcholslv(L,x);
    err_diff = norm(L'\(L\x)-y);
    assert(err_diff<1e-6, 'zchlslv_cpy is buggy');
end
elapsed_time = toc;

disp(['Test CHOLESKY    ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);