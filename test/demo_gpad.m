%% TEST GPAD IMPLEMENTATION:

clear all
clear classes
load ./archive/testCase1
NGpaditer = 5;
do_plots = (1==1);
do_monotone=true;
max_iter = 1000;
e_g = 1e-4; % infeasibility
e_V = 1e-4; % sub-optimality
N=10;


% SOLVE PROBLEM
% x+ = Ax + Bu +f
% xmin <= x <= xmax
% umin <= u <= umax
% cmin <= Fx + Gu <= cmax
% F=[]; G=[]; g=[]; % try with and without this line to make sure...
sys = LTISystem(A, B, [], F, G, [], g, xmin, xmax, umin, umax);

% 0.5(x'Qx+u'Ru+2u'S'x) + q'x + r'u 
% xN'*QN*xN + qN'*xN
mpc = MPCGPAD(sys);
mpc=mpc.setTerminalCost(dare(A, B, Q, R, []), []);
mpc.monotone=false;
mpc = mpc.setTerminalConstraints(FN,gN). ...
    setCost(Q, R).setHorizon(N);

[nx, nu, nc] = mpc.sys.getDimensions();
mpc.isInfHorizonCost = 1;
mpc = mpc.factor();
alpha = 1/mpc.calculateLipschitz();
cmap = jet(NGpaditer);

figure(1);
for j=1:NGpaditer
    % Test Dual Cost
    x0 = 1.2*rand(nx,1);    
    [y, diagnostics_simple]=mpc.solveDual(x0, alpha/1.001, max_iter, e_g, e_V);
    if (do_monotone==1),
        mpc.monotone=true;
        [y2, diagnostics_monotone]=mpc.solveDual(x0, alpha/1.001, max_iter, e_g, e_V);
        disp(['Termination : ' diagnostics_monotone.message ' | ' diagnostics_simple.message]);
    end
        
    if (do_plots)        
        col = cmap(j,:);
        subplot(311);
        if do_monotone==1, semilogy((max(diagnostics_monotone.primal_violation, e_g/10)),'--x','Color', col); hold on; end        
        semilogy(max(diagnostics_simple.primal_violation, e_g/10),'Color', col); hold on;
        grid on;
        if do_monotone==1, legend('Enforced Monotonicity','Normal Mode'); end
        title('Primal Infeasibility');
        
        subplot(312);hold on;
        grid on;
        title('Dual Cost');
        if do_monotone==1,  plot(diagnostics_monotone.dcosts/max(diagnostics_monotone.dcosts),'--x','Color',col); end;
        plot(diagnostics_simple.dcosts/max(diagnostics_monotone.dcosts),'Color',col);
        
        subplot(313);hold on;
        grid on;
        title('Duality Gap');
        if do_monotone==1, plot(diagnostics_monotone.dgap,'--x','Color',col); end
        plot(diagnostics_simple.dgap,'Color',col);

        if do_monotone==1,
            if ((diagnostics_monotone.iter > diagnostics_simple.iter))
                disp('Non-monotone wins');
            else
                disp('Monotone wins');
            end
        end
    end
end