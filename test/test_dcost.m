%TEST_DCOST Unit test function
%
% Tests the GPAD implementation of GPAD
%
%See also:
%dcost, pcost

%% FIRST PART - UNIT TEST: MPCGPAD#calculate_e(y,yN)
clear all
clear classes
load ./archive/randomPromlem1
tic;
Q = Q/1000;
A = 2*randn(nx,nx);
sys = LTISystem(A, B, [], F, G, [], g, xmin, xmax, umin, umax);
mpc = MPCGPAD(sys);
mpc = mpc.setTerminalConstraints(FN,gN).setCost(Q, R).setHorizon(8).factor();
[nx, nu, nc] = sys.getDimensions();
nf = size(mpc.FN,1);
np = nc + 2*nx + 2*nu;

for i=1:100
    if (i==10)
        load ./archive/testCase1
        sys = LTISystem(A, B, [], F, G, [], g, xmin, xmax, umin, umax);
        mpc = MPCGPAD(sys);
        mpc = mpc.setTerminalConstraints(FN,gN).setCost(Q, R).setHorizon(N).factor();
        [nx, nu, nc] = sys.getDimensions(); nf = size(mpc.FN,1); np = nc + 2*nx + 2*nu;
    end
    y = rand(np,mpc.N);
    yN = rand(nf,1);
    
    % Matrix `e` calculated by MPCGPAD#calculate_e(y,yN)
    e_mpc = mpc.calculate_e(y, yN);
    
    % Matrix `e` calculated from scratch:
    e_test = zeros(nx, mpc.N+1);
    e_test(:,mpc.N+1) = FN'*yN;
    Ctilde = [mpc.C eye(nx) -eye(nx) mpc.K' -mpc.K'];
    for k=mpc.N:-1:2
        e_test(:,k) = mpc.L*e_test(:,k+1) + Ctilde*y(:,k);
    end
    assert(norm(e_test-e_mpc,'inf')<1e-9, 'E is wrong');
end

elapsed_time = toc;
disp(['Test CALCUL. E   ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);


%% SECOND PART - TEST DCOST MEX (EXHAUSTIVE TEST)

clear all
tic;
nx=15;
nu=7;
gN = ones(nx,1);
N=8;
Q=0.001*eye(nx);
R=eye(nu);
for j=1:50
    A=rand(nx,nx); B=rand(nx,nu);
    FN = 10*rand*eye(nx);
    xmax=rand*ones(nx,1);     xmin = -xmax;
    umax=rand*ones(nu,1);     umin = -umax;
    nc=120;
    F=0.1*rand(nc,nx);
    G=rand(nc,nu);
    g=rand(nc,1);
    sys = LTISystem(A, B, [], F, G, [], g, xmin, xmax, umin, umax);
    mpc = MPCGPAD(sys);
    mpc = mpc.setTerminalConstraints(FN,gN).setCost(Q, R).setHorizon(N).factor();
    [nx, nu, nc] = sys.getDimensions();
    nf = size(mpc.FN,1);
    np = nc + 2*nx + 2*nu;
    
    x0=rand(nx,1);
    y = rand(np,mpc.N);
    yN = rand(nf,1);
    
    [Jdual_mpc, e_mpc, Useq_mpc, Xseq_mpc] = dcost(mpc, x0, y, yN);
    e_test = mpc.calculate_e(y,yN);
    % Check whether e_mpc was calculated correctly:
    assert(norm(e_mpc-e_test,'inf')<1e-10, 'Calculation of e_mpc failed');
    % Check whether Useq_mpc and Xseq_mpc are consistent:
    assert(norm(Xseq_mpc(:,1)-x0)<1e-10, 'Xseq(1)!=x0');
    for k=1:mpc.N
        assert(norm(A*Xseq_mpc(:,k)+B*Useq_mpc(:,k)-Xseq_mpc(:,k+1))<1e-10,...
            'Xseq and Useq are inconsistent');
    end
    % Calculate the dual cost:
    Jdual_test = 0;
    x=x0;
    Rbar = R + B'*mpc.QN*B;
    Rinv = inv(Rbar);
    Dtilde = [mpc.D zeros(nu,2*nx) -Rinv Rinv];
    
    for k=1:mpc.N
        u = mpc.K * x + Dtilde*y(:,k) + mpc.M * e_mpc(:,k+1);
        assert(norm(Useq_mpc(:,k)-u)<1e-7, ['Useq(' num2str(k) ') is wrong']);
        assert(norm(Xseq_mpc(:,k)-x)<1e-7, ['Xseq(' num2str(k) ') is wrong']);
        Jdual_test = Jdual_test + 0.5*(x'*Q*x + u'*R*u);
        gz = [F*x+G*u - g; x-xmax; xmin-x; u-umax; umin-u];
        Jdual_test = Jdual_test + y(:,k)' * gz;
        x = A*x + B*u;
    end
    Jdual_test = Jdual_test + 0.5*x'*mpc.QN*x + yN'*(FN*x-gN);
    assert(norm(x-Xseq_mpc(:,end))<1e-10, '');
    assert(abs(Jdual_mpc-Jdual_test)<1e-7, 'Dual Cost Discrepancy');
end
elapsed_time = toc;
disp(['Test DCOST#1     ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);

%% DETAILED TEST
clear all
load ./archive/testCase1
tic;
sys = LTISystem(A, B, [], F, G, [], g, xmin, xmax, umin, umax);

mpc = MPCGPAD(sys);
mpc = mpc.setTerminalConstraints(FN,gN).setCost(Q, R).setHorizon(8);
[nx, nu, nc] = mpc.sys.getDimensions();
mpc.isInfHorizonCost = 1;
mpc = mpc.factor();


% Test Dual Cost
x0 = rand(nx,1);
y = full(sprand(nc+2*(nx+nu),mpc.N,0.05));
yN = full(sprand(nf,1,0.05));


[dual_cost, eee, UU]=mpc.dcost(x0, y, yN);


C = mpc.C; K=mpc.K; L=mpc.L; s=mpc.s; D=mpc.D; M=mpc.M;
[A,B,f]=sys.getDynamics();
[F, G, cmin, cmax, umin, umax, xmin, xmax] = sys.getConstraints();
N=mpc.N;


Jdual = 0;
X=x0;
e = zeros(nx, N+1);
e(:, N+1) = mpc.FN'*yN;

Ctilde = [C eye(nx) -eye(nx) K' -K'];
CtildeY = zeros(nx, 1);
for k=N:-1:2,
    e_test = L * e(:, k+1) +  Ctilde*y(:,k);
    e(:, k) = L * e(:, k+1);
    i=0;
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
    if ~isempty(umax),
        CtildeY=CtildeY-K'*y(i+1:i+nu,k);
        doAddCtildeY = true;
    end
    if (doAddCtildeY), e(:, k) = e(:, k) + CtildeY; end;
    assert(norm(e_test - e(:,k))<1e-10);
end % end of the first loop
assert(norm(eee-e)<1e-10);
e_mpc = mpc.calculate_e(y,yN);
assert(norm(eee-e_mpc)<1e-10);
assert(norm(e-e_mpc)<1e-10);

Rbar = mpc.R+B'*mpc.QN*B;
Rinv = inv(Rbar);
Dtilde=[-Rinv*G' zeros(nu, 2*nx) -Rinv Rinv];
for k=1:N,
    U = K*X  + M * e(:,k+1); i=0;
    Utest = U + Dtilde*y(:,k);
    if ~isempty(cmax), U = U + D*y(1:nc,k); i=i+nc;  end
    if ~isempty(cmin), U = U - D*y(i+1:i+nc,k); i=i+nc;  end
    i = i + nx*(~isempty(xmin)+~isempty(xmax));
    if ~isempty(umax), U = U - mpc.RbarChol'\(mpc.RbarChol\y(i+1:i+nu,k)); i=i+nu;  end
    if ~isempty(umin), U = U + mpc.RbarChol'\(mpc.RbarChol\y(i+1:i+nu,k)); end
    Jdual = Jdual + 0.5*(X'*mpc.Q*X + U'*mpc.R*U);
    Jdual = Jdual + y(:,k)' * primalFeasibility(mpc, X, U);
    X = A*X + B*U;
    assert(norm(U-Utest)<1e-10);
    assert(norm(Utest-UU(:,k))<1e-10);
end
Jdual = Jdual + 0.5*X'*mpc.QN*X;
Jdual = Jdual + yN'*(mpc.FN*X-mpc.gN);
error = Jdual-dual_cost;
assert(abs(error)<1e-10);
elapsed_time = toc;

disp(['Test DUAL COST   ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);


%% TEST - Addume F, G, cmin, cmax are empty
clear all
tic;
nx=15;
nu=7;
gN = ones(nx,1);
N=8;
Q=0.001*eye(nx);
R=eye(nu);

A=rand(nx,nx); B=rand(nx,nu); FN = 10*rand*eye(nx);
xmax=rand*ones(nx,1);xmin = -xmax; umax=rand*ones(nu,1);umin = -umax;
sys = LTISystem(A, B, [], [], [], [], [], xmin, xmax, umin, umax);
mpc = MPCGPAD(sys);
mpc = mpc.setTerminalConstraints(FN,gN).setCost(Q, R).setHorizon(N).factor();
[nx, nu, nc] = sys.getDimensions();
nf = size(mpc.FN,1);
np = nc + 2*nx + 2*nu;

x0=rand(nx,1);
y = rand(np,mpc.N);
yN = rand(nf,1);


[d, Ustar, Xstar] = dcost(A, B, [], Q, R, [], [], [], mpc.QN, [], FN, gN, ...
    [], [], [], [], xmin, xmax, umin, umax, mpc.D, mpc.L, mpc.M, mpc.K, ...
    mpc.C, mpc.s, mpc.d,  mpc.RbarChol, y, yN, x0, N);


[Jdual_mpc, e_mpc, Useq_mpc, Xseq_mpc] = dcost(mpc, x0, y, yN);

assert(abs(Jdual_mpc-d)<1e-7, 'Jdual not calculated propertly in the case cmin=cmax=[]');

elapsed_time = toc;
disp(['Test DCOST_EMPTY ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);
