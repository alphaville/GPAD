function [K, L, D, M, C, s, d, e, Jdual, Xstar, Ustar] = ...
    solve_dual(A, B, f, Q, R, S, q, r, QN, qN, F, G, cmin, cmax, xmin, xmax, umin, umax, FN, gN, N, x0, y, yN)

nx=size(A,1);
nu=size(B,2);
nc=size(F,1);

%% FACTOR STEP
P = QN;
BtransP = B'*P;
Rbar = R  +  BtransP * B;
Sbar = BtransP * A;
diff_P=norm(Q+A'*P*A-(Sbar'*(Rbar\Sbar))-P);
assert(diff_P<1e-6, ...
    ['Problem with QN - Difference = ' num2str(diff_P)]);
if ~isempty(S), Sbar = Sbar + S; end
RbarChol = chol(Rbar,'lower');
K = -RbarChol'\(RbarChol\Sbar);
D = []; if ~isempty(G), D = -RbarChol'\(RbarChol\G'); end
M = -RbarChol'\(RbarChol\B');
if (isempty(r) && isempty(f)), d=[]; else
    dtemp = zeros(nx,1);
    if ~isempty(f), dtemp = BtransP * f; end
    if ~isempty(r), dtemp = dtemp + r; end
    d = -RbarChol'\(RbarChol\dtemp);
end
L = (A + B*K)';
if ~isempty(F), C = (F + G * K)'; end
s=[];
if ~isempty(q) || ~isempty(r) || ~isempty(f),
    s = zeros(nx,1);
    if ~isempty(r), s = s + K * r; end
    if ~isempty(f), s = s + L * P * f; end
    if ~isempty(q), s = s + q; end
end

%% CALCULATE e
e = zeros(nx, N+1);
e(:,N+1) = FN' * yN;
if ~isempty(qN), e(:,N+1) = e(:,N+1) + qN; end
for k=N:-1:2,
    e(:, k) = L * e(:, k+1);
    if ~isempty(s), e(:, k) = e(:, k) + s; end
    i=0;
    CtildeY = zeros(nx, 1);
    doAddCtildeY = false;
    if ~isempty(cmax),
        CtildeY=C*y(1:nc,k); i=i+nc; doAddCtildeY = true;
    end
    if ~isempty(cmin),
        CtildeY=CtildeY-C*y(i+1:i+nc,k);
        i=i+nc; doAddCtildeY = true;
    end
    if ~isempty(xmax),
        CtildeY=CtildeY+y(i+1:i+nx,k);
        i=i+nx; doAddCtildeY = true;
    end
    if ~isempty(xmin),
        CtildeY=CtildeY-y(i+1:i+nx,k);
        i=i+nx; doAddCtildeY = true;
    end
    if ~isempty(umax),
        CtildeY=CtildeY+K'*y(i+1:i+nu,k);
        i=i+nu; doAddCtildeY = true;
    end
    if ~isempty(umin),
        CtildeY=CtildeY-K'*y(i+1:i+nu,k);
        doAddCtildeY = true;
    end
    if (doAddCtildeY), e(:, k) = e(:, k) + CtildeY; end;
end

%% SOLUTION STEP
Jdual = 0;
X=x0;
Ustar=zeros(nu, N);
Xstar=zeros(nx, N+1);
Xstar(:,1)=x0;
for k=1:N,
    U = K*X  + M * e(:,k+1); i=0;     
    if ~isempty(cmax), U = U + D*y(1:nc,k); i=i+nc;  end
    if ~isempty(cmin), U = U - D*y(i+1:i+nc,k); i=i+nc;  end
    i = i + nx*(~isempty(xmin)+~isempty(xmax));
    if ~isempty(umax), U = U - RbarChol'\(RbarChol\y(i+1:i+nu,k)); i=i+nu;  end
    if ~isempty(umin), U = U + RbarChol'\(RbarChol\y(i+1:i+nu,k)); end
    Jdual = Jdual + 0.5*(X'*Q*X + U'*R*U);
    if ~isempty(S), Jdual = Jdual + U'*S*X; end
    if ~isempty(cmax) || ~isempty(cmin), FxGu = F*X + G*U; end
    np = nc*(~isempty(cmax) + ~isempty(cmin)) + ...
        nx*(~isempty(xmax)+ ~isempty(xmin)) + ...
        nu*(~isempty(umax)+ ~isempty(umin));
    g = zeros(np,1); i=0;
    if ~isempty(cmax), g(1:nc) = FxGu - cmax; i = i+nc; end
    if ~isempty(cmin), g(i+1:i+nc) = cmin - FxGu; i = i+nc; end
    if ~isempty(xmax), g(i+1:i+nx) = X - xmax; i = i+nx; end
    if ~isempty(xmin), g(i+1:i+nx) = xmin - X; i = i+nx; end
    if ~isempty(umax), g(i+1:i+nu) = U - umax; i = i+nu; end
    if ~isempty(umin), g(i+1:i+nu) = umin - U; end
    Jdual = Jdual + y(:,k)' * g;
    X = state_update(A, B, f, X, U);
    Ustar(:,k) = U;
    Xstar(:,k+1) = X;
end
Jdual = Jdual + 0.5*X'*QN*X;
Jdual = Jdual + yN'*(FN*X-gN);

end