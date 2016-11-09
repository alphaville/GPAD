%DEMO_MPC demonstrates the use of the GPADToolbox and the classes LTISystem
%and MPCGPAD to design an MPC feedback controller.
%
%See also:
%demo_gpad
clear all
clear classes
load ./archive/testCase2


% GPAD Data:
max_iter = 6000;
e_g = 2e-4; % infeasibility
e_V = 1e-4; % sub-optimality
N=40; % Prediction Horizon
R=5*R;
% Initial state:
x0 = rand(nx,1);

% Define the LTI system:
sys = LTISystem(A, B, [], [], [], [], [], xmin, xmax, umin, umax);

% Construct the controller:
mpc = MPCGPAD(sys); % Define a controller for the above LTI system
Pf=dare(A, B, Q, R, []); % Calculate the terminal cost
mpc=mpc.setTerminalCost(Pf,[]); % Set the terminal cost Vf(x) = x' Pf x
mpc.monotone=true; % Enforce monotonicity (true/false)
mpc = mpc.setTerminalConstraints(FN,gN). ...
    setCost(Q, R).setHorizon(N); % Set the terminal constraints

[nx, nu, nc] = mpc.sys.getDimensions();
mpc.isInfHorizonCost = 1; % The specified cost is the infinite horizon cost
mpc = mpc.factor(); % Perform the factorization (Important step)
alpha = 1/mpc.calculateLipschitz(); % Calculate a Lipschitz constant
%% Simulate
Nsim=120;
X=zeros(nx, Nsim+1);
U=zeros(nu, Nsim);
X(:,1)=x0;
sys=sys.setState(x0);
Diagnostics = struct();
for k=1:Nsim,
    [Ustar, Xstar, Diag] = mpc.control(sys.getState, alpha, e_g, e_V, max_iter);
    sys = sys.stateupdate(Ustar(:,1));    
    U(:, k)=Ustar(:,1);
    X(:,k+1)=sys.getState();    
    Diagnostics(k).Diag=Diag;
end

%% Plot
figure
hold on;
subplot(211);
plot(X'); hold on; ylabel('State'); grid on
subplot(212);
plot(U'); hold on; ylabel('Input');grid on
