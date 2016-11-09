clear all
load ./archive/testCase1.mat
tic;
Q=rand(nx,nx); Q=Q'*Q;
Q=Q/1000;
R=rand(nu,nu);R=R'*R;
Pf=dare(A, B, Q, R, []);
x0=rand(nx,1);
np = size(F,1) + 2*nx + 2*nu;
w=rand(np,N);
wN=rand(nf,1);
[K, L, D, M, C, s, d, e, Jdual, Xstar, Ustar] = ...
    solve_dual(A, B, [], Q, R, [], [], [], Pf, [], ...
    F, G, [], g, xmin, xmax, umin, umax, FN, gN, N, x0, w, wN);

Rbar = R  +  B'*Pf*B;
Rchol = chol(Rbar,'lower');

slack_test = zeros(np, N);
e_test = calculate_e(FN, [], L, C, K, s, 1, 0, 0, 0, 0, 0, N, w, wN);
U_test = zeros(nu, N);
X_test = zeros(nx, N+1);
X_test(:,1)=x0;
for k=1:N,
    U_test(:,k) = K*X_test(:,k)  + M * e_test(:,k+1); i=0;
    U_test(:,k) = U_test(:,k) + D*w(1:nc,k); i=i+nc;
    i = i + nx*2;
    U_test(:,k) = U_test(:,k) - Rchol'\(Rchol\w(i+1:i+nu,k)); i=i+nu; 
    U_test(:,k) = U_test(:,k) + Rchol'\(Rchol\w(i+1:i+nu,k)); 
    slack_test(:,k) = [F*X_test(:,k)+G*U_test(:,k)-g;
        X_test(:,k)-xmax;
        xmin-X_test(:,k);
        U_test(:,k)-umax;
        umin-U_test(:,k)];
    X_test(:,k+1) = state_update(A,B,[],X_test(:,k), U_test(:,k));
end
slackN_test = FN*X_test(:,N+1)-gN;

[slack, slackN, Xstar, Ustar, y, yN]=dgrad(A, B, [], Q, R, [], [], [], Pf, [], FN, gN, F, G, [], g, xmin, ...
    xmax, umin, umax, D, L, M, K, C, s, d, Rchol, w, wN, x0, N, alpha);

assert(norm(slack-slack_test)<1e-7, 'slack not calculated properly');
assert(norm(slackN-slackN_test)<1e-7, 'slackN not calculated properly');
assert(norm(Xstar-X_test)<1e-7, 'Xstar not calculated properly');
assert(norm(Ustar-U_test)<1e-7, 'Ustar not calculated properly');

elapsed_time = toc;
disp(['Test D. GRADIENT ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);
