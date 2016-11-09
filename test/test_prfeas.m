clear all
clear classes
load ./archive/testCase1

tic;
for i=1:100
    x=rand(nx,1);
    u=rand(nu,1);
    
    % SOLVE PROBLEM
    % x+ = Ax + Bu +f
    % xmin <= x <= xmax
    % umin <= u <= umax
    % cmin <= Fx + Gu <= cmax
    
    mpc = MPCGPAD(LTISystem(A, B, [], F, G, -g, g, xmin, xmax, umin, umax));
    gxu = mpc.primalFeasibility(x,u);
    gxu_gpad = primal_feasibility(F,G,-g,g,umin,umax,xmin,xmax,x,u);
    assert(norm(gxu-gxu_gpad)<1e-6, 'primal_feasibility: Wrong result');
    
    mpc = MPCGPAD(LTISystem(A, B, [], F, G, [], g, [], [], [], []));
    gxu = mpc.primalFeasibility(x,u);
    gxu_gpad = primal_feasibility(F,G,[],g,[],[],[],[],x,u);
    assert(norm(gxu-gxu_gpad)<1e-6, 'primal_feasibility: Wrong result (2)');
    
    gxu_gpad = primal_feasibility([],[],[],[],[],[],[],[],x,u);
    assert(isempty(gxu_gpad), 'gxu_gpad should be empty if no constraints are given');
end
elapsed_time=toc;
disp(['Test PR. FEASIB. ~ OK.  Completed in ' ...
    num2str(round(100*elapsed_time)/100) 's']);