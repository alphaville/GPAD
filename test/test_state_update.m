NStateUpdateIter=1000;
tic;
for i=1:NStateUpdateIter
    nx=1+floor(120*rand);
    nu=1+floor(80*rand);
    A=rand(nx,nx); 
    B=rand(nx,nu); 
    f=[]; 
    if nx>10, f=rand(nx,1); end
    x=rand(nx,1); 
    u=rand(nu,1);
    xnew = state_update(A,B,f,x,u);
    xnew_ = A*x + B*u;
    if ~isempty(f), xnew_ = xnew_ + f; end
    assert(norm(xnew-xnew_)<1e-10);
end
elapsed_time=toc;
disp(['Test STATEUPDATE ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);