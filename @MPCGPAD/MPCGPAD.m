classdef MPCGPAD
    %MPCGPAD is an implementation of the Accelerated Dual
    %Gradient-Projection Algorithm for Linear Model Predictive Control of
    %LTI and LTV systems.
    %
    %Part of the GPAD-Toolbox.
    %
    %References:
    % [1] P. Patrinos and A. Bemporad, "An Accelerated Dual
    % Gradient-Projection Algorithm for Embedded Linear Model Predictive
    % Control", 2012, Under review.
    
    properties
        % The dynamical system which the MPC controller needs to stabilise.
        sys=DynamicalSystem();
        
        % Identify the structure of the problem:
        %  (TODO: Maybe we need to make all these fields private)
        
        % Whether the terminal cost is equal to the infinite horizon cost
        % which is calculated by a Ricatti iteration.
        isInfHorizonCost = false;
    end
    
    properties
        %(Access=private)
        % N is the prediction horizon of the system.
        N=[];
        monotone=true;
       
        % Matrices of the stage cost function:
        % l(x,u) = 1/2 z' [Q S'; S R] z + [q;r]' z, where z=[x;u].
        % R is symmetric and positive definite
        % [Q S'; S R] is symmetric and positive definite
        Q=[]; S=[]; R=[]; q=[]; r=[];
        
        % Matrices of the terminal cost: Vf(x) = 1/2 x' QN x + qN' x
        % QN is symmetric and positive semidefinite.
        QN=[]; qN=[];
        
        % Terminal constraints: FN*xN <= gN
        FN=[]; gN=[];
        
        % Matrices of the factorisation:
        K=[]; D=[]; M=[]; d=[]; L=[];
        C=[]; s=[]; y=[];
        RbarChol=[];
        nf=[];
        
        % 0: dense, 1: Q=aI, 2: Q=diag, 3: Q=sparse
        classQ=0; classR=0;
        
        is_q_zero=true;
        is_r_zero=true;
        is_S_zero=true;
    end % end of private properties
    
    methods
        function obj = MPCGPAD(varargin)
            %MPCGPAD is the constructor method for the class MPCGPAD.
            %Objects can be constructed in the following ways:
            %
            % - mpc = MPCGPAD();
            %  This will create an empty MPCGPAD object. The parameters of
            %  the MPC problem can be specified later.
            %
            % - mpc = MPCGPAD(sys);
            %  Where sys is an LTI or and LTV system. In that case, the
            %  cost matrices are set to the default values Q=eye(nx) and
            %  R=eye(nu) and QN is set to the be the optimal cost of the
            %  infinite-horizon problem. This is determined as the solution
            %  of the following Ricatti-type equation:
            %   P = Q + A'PA - A'PB(R + B'PB)B'PA
            %
            if (nargin==1 && isa(varargin{1},'LTISystem'))
                obj.sys = varargin{1};
                obj.N=10;
                [nx, nu] = obj.sys.getDimensions();
                obj.Q = eye(nx); obj.classQ = 1;
                obj.R = eye(nu); obj.classR = 1;
                obj.is_q_zero = true;
                obj.is_r_zero = true;
                obj.isInfHorizonCost=true;
                [A, B]=obj.sys.getDynamics();                 
                obj.QN=dare(A, B, obj.Q, obj.R, obj.S);
                
                % diff_QN = norm(obj.Q+A'*obj.QN*A-((B'*obj.QN*A)'*((obj.R+B'*obj.QN*B)\(B'*obj.QN*A)))-obj.QN);
                % assert(norm(diff_QN)<1e-5);
            end
            if isempty(obj.S)
                obj.is_S_zero = true;
            end
        end
        
        function obj = setHorizon(obj, N)
            %SETHORIZON is a setter method for the prediction horizon.
            obj.N=N;
        end
        
        function N = getHorizon(obj)
            %GETHORIZON returns the prediction horizon.
            N = obj.N;
        end
        
        function [Q, R, S, q, r] = getCost(obj)
            %GETCOST returns the cost matrices of the problem.
            Q=obj.Q;
            R=obj.R;
            S=obj.S;
            q=obj.q;
            r=obj.r;
        end
        
        function obj = setCost(obj, Q, R, S, q, r)
            %SETCOST is a setter method for the matrices of the cost
            %function of the problem.
            %
            %See also:
            %MPCGPAD/setTerminalCost
            obj.Q=Q;
            obj.R=R;
            if nargin>=4, obj.S=S; end
            if nargin>=5, obj.q=q; end
            if nargin>=6, obj.r=r; end
        end
        
        function obj = setTerminalCost(obj, QN, qN)
            %SETTERMINALCOST is a setter method for the terminal cost
            %matrices of the MPC problem.
            obj.QN=QN;
            obj.qN=qN;
        end
        
        function [QN, qN] = getTerminalCost(obj)
            %GETTERMINALCOST returns the terminal cost matrices.
            QN = obj.QN;
            qN = obj.qN;
        end
        
        function obj = setTerminalConstraints(obj, FN, gN)
            %SETTERMINALCONSTRAINTS is a setter method for the terminal
            %cosntraints' matrices FN and gN, such that FN*xN <= gN.
            obj.FN=FN;
            obj.gN=gN;
            obj.nf = size(FN,1);
            assert(obj.nf==size(gN,1),'FN and gN have incompatible sizes');
        end
        
        function [FN, gN] = getTerminalConstraints(obj)
            %GETTerminalConstraints returns the terminal constraints'
            %matrices FN and gN such that xN satisfies: FX*xN<=gN.
            FN = obj.FN;
            gN = obj.gN;
        end
        
        function L = calculateLipschitz(obj)
            [F, G] = obj.sys.getConstraints();
            [nx, nu] = obj.sys.getDimensions();
            FGtilde = [ F   G
                       eye(nx) zeros(nx,nu)
                        -eye(nx) zeros(nx,nu)
                        zeros(nu,nx) eye(nu)
                        zeros(nu,nx) -eye(nu)];
                    if (~isempty(obj.S))
                        Qtilde = [obj.Q obj.S'
                                  obj.S obj.R];
                    else
                        Qtilde = blkdiag(obj.Q, obj.R);
                    end
                   L = norm(FGtilde)^2/min(min(eig(Qtilde)), min(eig(obj.QN)));
        end % end of method: calculateLipschitz
        
        function J = pcost(obj, X, U)
            J = pcost(obj.Q, obj.R, obj.S, obj.q, obj.r, ...
                obj.QN, obj.qN, obj.N, X, U);
        end
        
        function e = calculate_e(obj, w, wN)
            [nx, nu, nc] = obj.sys.getDimensions();
            [~, ~, cmin, cmax, umin, umax, xmin, xmax] = obj.sys.getConstraints();
            e = zeros(nx, obj.N+1);
            e(:,obj.N+1) = obj.FN' * wN;
            if ~isempty(obj.qN), e(:,obj.N+1) = e(:,obj.N+1) + obj.qN; end
            for k=obj.N:-1:2,
                e(:, k) = obj.L * e(:, k+1);
                if ~isempty(obj.s), e(:, k) = e(:, k) + obj.s; end
                i=0;
                CtildeY = zeros(nx, 1);
                doAddCtildeY = false;
                if ~isempty(cmax),
                    CtildeY=obj.C*w(1:nc,k); i=i+nc; doAddCtildeY = true;
                end
                if ~isempty(cmin),
                    CtildeY=CtildeY-obj.C*w(i+1:i+nc,k);
                    i=i+nc; doAddCtildeY = true;
                end
                if ~isempty(xmax),
                    CtildeY=CtildeY+w(i+1:i+nx,k);
                    i=i+nx; doAddCtildeY = true;
                end
                if ~isempty(xmin),
                    CtildeY=CtildeY-w(i+1:i+nx,k);
                    i=i+nx; doAddCtildeY = true;
                end
                if ~isempty(umax),
                    CtildeY=CtildeY+obj.K'*w(i+1:i+nu,k);
                    i=i+nu; doAddCtildeY = true;
                end
                if ~isempty(umin),
                    CtildeY=CtildeY-obj.K'*w(i+1:i+nu,k);
                    doAddCtildeY = true;
                end
                if (doAddCtildeY), e(:, k) = e(:, k) + CtildeY; end;
            end % end first for loop - e(2:end, k)
        end % end of method: calculate_e 
        
        function [Ustar, Xstar, Diag] = control(obj, x0, alpha, e_g, e_V, max_iter)
            %CONTROL computes the control action given the state of the
            %system, the maximum primal constraints' violation and the
            %maximum suboptimality.
            %
            %Syntax:
            % [Ustar, Xstar, Diag] = mpc.CONTROL(x0, alpha, e_g, e_V, max_iter)
            %
            %Input Arguments:
            %~ mpc      : An MPCGPAD object
            %~ x0       : The current state
            %~ alpha    : The reciprocal of a Lipschitz constant for the
            %             dual gradient which can be calculated using the
            %             method calculateLipschitz
            %~ e_g      : Maximum permissible constraint violation
            %~ e_V      : Primal suboptimality bound
            %~ max_iter : Maximum Iterations
            %
            %Output Arguments:
            %~ Ustar    : The (e_g, e_V)-optimal sequence of inputs
            %~ Xstar    : The corresponding sequence of states
            %~ Diag     : Diagnostics
            %
            %See also:
            % MPCGPAD/calculateLipschitz
            
            [A, B, f] = obj.sys.getDynamics();
            [F, G, cmin, cmax, umin, umax, xmin, xmax] = obj.sys.getConstraints();
            [Ustar, Xstar, Diag] = gpad(A, B, f, obj.Q, obj.R, obj.S, obj.q, obj.r,...
                obj.QN, obj.qN, obj.FN, obj.gN, F, G, cmin, cmax, xmin, xmax, umin, ...
                umax, obj.D, obj.L, obj.M, obj.K, obj.C, obj.s, obj.d, ...
                obj.RbarChol, x0, obj.N, alpha, e_g, e_V, max_iter);
        end
        
        function [y, diagnostics] = ...
                solveDual(obj, x0, alpha, maxiter, epsilon_g, epsilon_V)
            
            % System data:
            [F, G, cmin, cmax, umin, umax, xmin, xmax] = obj.sys.getConstraints();
            [A, B, f] = obj.sys.getDynamics();
            [nx, nu, nc] = obj.sys.getDimensions();
            np = nc*(~isempty(cmax) + ~isempty(cmin)) + ...
                nx*(~isempty(xmax)+ ~isempty(xmin)) + ...
                nu*(~isempty(umax)+ ~isempty(umin));
                        
            % Initial Conditions:
            check2=0;check3=0;check4=0;
            if (nargout>=2), 
                dcosts = zeros(maxiter,1); diagnostics.message='max iter reached'; 
            end
            loop=true; iter=0; theta=1; theta_p=1;
            y = zeros(np,obj.N);
            yN = zeros(obj.nf,1);     
            y_p=y; yN_p=yN;
                
            slack_bar = zeros(np, obj.N);
            slackN_bar = zeros(obj.nf, 1);
            
            if (nargout>=2),
                diagnostics.primal_violation = zeros(maxiter,1);
                diagnostics.dgap = zeros(maxiter,1);
            end
            
            %% THE LOOP: 
            while loop && iter<maxiter
                % Extrapolation:
                beta = theta*(1/theta_p-1);
                
                w = y + beta*(y-y_p);
                wN = yN + beta*(yN-yN_p);
                               
                y_p=y;
                yN_p=yN;
                                 
                %% Calculate the Gradient:                                
                [slack, slackN, X, U, y, yN]=dgrad(A, B, f, obj.Q, obj.R, ...
                    [], [], [], obj.QN, [], obj.FN, obj.gN, F, G, cmin, cmax, xmin, ...
                    xmax, umin, umax, obj.D, obj.L, obj.M, obj.K, ...
                    obj.C, obj.s, obj.d, obj.RbarChol, w, wN, x0, obj.N, alpha);
                                              
                
                %% Check if the gradient was calculated correctly:                                                                                                                          
                % Monotonicity Enforcement:                   
                if (obj.monotone && iter>1)
                    if (obj.dcost(x0,y_p,yN_p)>obj.dcost(x0,y,yN)),
                        y=y_p; yN=yN_p;
                    end
                end
                
                max_viol = max(max(max(slack)), max(slackN));
                if (nargout>=2), diagnostics.primal_violation(iter+1)=max_viol; end
                
                % Calculate slack_bar
                slack_bar  = (1-theta)*slack_bar+theta*slack;
                slackN_bar = (1-theta)*slackN_bar+theta*slackN;
                
                max_viol_bar = max(max(max(slack_bar)),max(slackN_bar));
                if (nargout>=2), dcosts(iter+1) = obj.dcost(x0, y, yN); end
                
                %% TERMINATION CRITERION:
                if (nargout==2)
                    for k=1:obj.N, 
                        diagnostics.dgap(iter+1) = diagnostics.dgap(iter+1) - ...
                            w(:,k)'*slack(:,k); 
                    end
                end
                diagnostics.dgap(iter+1) = diagnostics.dgap(iter+1)-wN'*slackN;                        
                if max_viol_bar<=epsilon_g, % POINT #1
                    diagnostics.message='max_viol_bar<eg'; 
                    loop = false;
                elseif max_viol<=epsilon_g % POINT %2
                    if min(wN)>=0  && min(min(w))>=0
                        check2=check2+1; gap = 0;
                        for k=1:obj.N, gap = gap - w(:,k)'*slack(:,k); end
                        gap = gap-wN'*slackN;
                        if gap<=epsilon_V, diagnostics.message='gap<eV'; % POINT #3
                            loop = false; 
                        else % POINT 4
                            check3=check3+1; primal_cost = obj.pcost(X,U);
                            if gap<=epsilon_V/(1+epsilon_V)*primal_cost
                                diagnostics.message='gap<eV/(1+eV)*primal_cost'; 
                                loop = false;
                            end
                        end
                    else% if min(wN)>=0  && min(min(w))>=0
                        check4=check4+1; primal_cost = obj.pcost(X,U);
                        dual_cost = obj.dcost(x0, y, yN);
                        rel_gap = (primal_cost-dual_cost)/max(dual_cost,1);
                        if rel_gap<=epsilon_V, 
                            diagnostics.message='rel_gap<eV'; loop = false; 
                        end
                    end                    
                end % END: Termination Criterion
                
                
                %% UPDATE :                
                %Update Theta:
                theta_p = theta;
                theta2 = theta^2;
                theta = (sqrt(theta2^2+4*theta2)-theta2)/2;
                iter = iter + 1;
            end % end of major while loop
            if (nargout>=2)
                diagnostics.dcosts=dcosts(1:iter);
                diagnostics.dgap=diagnostics.dgap(1:iter);
                diagnostics.primal_violation=diagnostics.primal_violation(1:iter);
                diagnostics.check=[check2;check3;check4];   
                diagnostics.iter=iter;
            end
        end % end of method: solveDual                
        
        function [Jdual, e, UU, XX] = dcost(obj, x0, y, yN)
            %DCOST returns the optimal dual cost. It is an implementation
            %of Algorithm 4 in [1].
            Jdual = 0; X=x0;
            [nx, nu, nc] = obj.sys.getDimensions();            
            [A, B, f] = obj.sys.getDynamics(); UU=zeros(nu, obj.N);
            XX=zeros(nx, obj.N+1);
            XX(:,1)=x0;
            [~, ~, cmin, cmax, umin, umax, xmin, xmax] = obj.sys.getConstraints();
            % First for loop - Calculate e(k) for k=2,...,N           
            e = obj.calculate_e(y, yN);
            for k=1:obj.N,
                U = obj.K*X  + obj.M * e(:,k+1); i=0;
                if ~isempty(cmax), U = U + obj.D*y(1:nc,k); i=i+nc;  end
                if ~isempty(cmin), U = U - obj.D*y(i+1:i+nc,k); i=i+nc;  end
                i = i + nx*(~isempty(xmin)+~isempty(xmax));
                if ~isempty(umax), U = U - obj.RbarChol'\(obj.RbarChol\y(i+1:i+nu,k)); i=i+nu;  end
                if ~isempty(umin), U = U + obj.RbarChol'\(obj.RbarChol\y(i+1:i+nu,k)); end
                Jdual = Jdual + 0.5*(X'*obj.Q*X + U'*obj.R*U);
                if ~isempty(obj.S), Jdual = Jdual + U'*obj.S*X; end
                Jdual = Jdual + y(:,k)' * obj.primalFeasibility(X,U);
                X = state_update(A, B, f, X, U);
                UU(:,k)=U;
                XX(:,k+1)=X;
            end
            Jdual = Jdual + 0.5*X'*obj.QN*X;
            Jdual = Jdual + yN'*(obj.FN*X-obj.gN);
        end % end of method: dcost
        
        function g = primalFeasibility(obj, x, u)
            [F, G, cmin, cmax, umin, umax, xmin, xmax] = obj.sys.getConstraints();
            [nx, nu, nc] = obj.sys.getDimensions();
            if ~isempty(cmax) || ~isempty(cmin), FxGu = F*x + G*u; end
            np = nc*(~isempty(cmax) + ~isempty(cmin)) + ...
                nx*(~isempty(xmax)+ ~isempty(xmin)) + ...
                nu*(~isempty(umax)+ ~isempty(umin));
            g = zeros(np,1); i=0;
            if ~isempty(cmax), g(1:nc) = FxGu - cmax; i = i+nc; end
            if ~isempty(cmin), g(i+1:i+nc) = cmin - FxGu; i = i+nc; end
            if ~isempty(xmax), g(i+1:i+nx) = x - xmax; i = i+nx; end
            if ~isempty(xmin), g(i+1:i+nx) = xmin - x; i = i+nx; end
            if ~isempty(umax), g(i+1:i+nu) = u - umax; i = i+nu; end
            if ~isempty(umin), g(i+1:i+nu) = umin - u; end
        end % end of method: primalFeasibility
        
        function obj = factor(obj)
            %FACTOR is an implementation of Algorithm 3 in [1].
            %It returns the matrices K, D, M, d, L, C and s of Algorithm 3.
            %Unless QN (the terminal cost matrix) is the infinite horizon
            %cost, these matrices are 3D matrices so that K_k = K(:,:,k)
            %for k=1,...,N.
            [A, B, f] = obj.sys.getDynamics();
            [F, G] = obj.sys.getConstraints();
            [nx, nu, nc] = obj.sys.getDimensions();
            isf0 = obj.sys.isfzero();
            if obj.isInfHorizonCost,
                P = obj.QN;
                BtransP = B'*P;
                Rbar = obj.R  +  BtransP * B;
                Sbar = BtransP * A;
                if ~obj.is_S_zero,
                    Sbar = Sbar + obj.S;
                end
                % Rbar is positive definite and symmetric - we'll use the
                % Cholesky factorisation.
                obj.RbarChol = chol(Rbar,'lower');
                obj.K = -obj.RbarChol'\(obj.RbarChol\Sbar);
                obj.D = [];
                if ~isempty(G),
                    obj.D = -obj.RbarChol'\(obj.RbarChol\G');
                end
                obj.M = -obj.RbarChol'\(obj.RbarChol\B');
                if (obj.is_r_zero && isf0)
                    obj.d=[];
                else
                    dtemp = zeros(nx,1);
                    if ~obj.sys.isfzero(),
                        dtemp = BtransP * f;
                    end
                    if ~obj.is_r_zero,
                        dtemp = dtemp + obj.r;
                    end
                    obj.d = -obj.RbarChol'\(obj.RbarChol\dtemp);
                end
                obj.L = (A + B*obj.K)';
                if ~isempty(F)
                    obj.C = (F + G * obj.K)';
                end
                if ~obj.is_q_zero || ~obj.is_r_zero || ~isf0,
                    obj.s = zeros(nx,1);
                    if ~obj.is_r_zero, obj.s = obj.s + obj.K * obj.r; end
                    if ~isf0, obj.s = obj.s + obj.L * P * f; end
                    if ~obj.is_q_zero, obj.s = obj.s + obj.q; end
                else
                    obj.s=[];
                end
            else % The terminal cost is not the infinite horizon cost:
                P=obj.QN;
                obj.M = zeros(nu, nx, obj.N);
                obj.K = zeros(nu, nx, obj.N);
                obj.D = zeros(nu, nc, obj.N);
                obj.L = zeros(nx, nx, obj.N);
                obj.C = zeros(nx, nc, obj.N);
                if ~obj.is_q_zero || ~obj.is_r_zero || ~isf0,
                    obj.s = zeros(nx, 1, obj.N);
                else
                    obj.s=[];
                end
                if ~obj.is_r_zero || ~isf0,
                    obj.d = zeros(nu, 1, obj.N);
                else
                    obj.d=[];
                end
                for k=obj.N:-1:1
                    BtP  = B'*P;
                    Rbar = obj.R + BtP*B;
                    Sbar = BtP*A;
                    obj.RbarChol = chol(Rbar, 'lower');
                    obj.K(:,:,k)  = -obj.RbarChol'\(obj.RbarChol\Sbar);
                    obj.D(:,:,k)  = -obj.RbarChol'\(obj.RbarChol\G');
                    obj.M(:,:,k)  = -obj.RbarChol'\(obj.RbarChol\B');
                    obj.L(:,:,k)  = (A + B*obj.K(:,:,k))';
                    obj.C(:,:,k)  = (F + G*obj.K(:,:,k))';
                    P         = obj.Q + A'*P*A + Sbar'*obj.K(:,:,k);
                    if ~obj.is_q_zero, obj.s(:,1,k) = obj.q; end
                    if ~obj.is_r_zero,
                        obj.s(:,1,k)=obj.K(:,:,k)*obj.r;
                        obj.d = obj.r;
                    end
                    if ~obj.sys.is_f_zero,
                        obj.s(:,1,k)=obj.L(:,:,k)*P*f;
                        obj.d = obj.d + BtP * f;
                    end
                    if ~isempty(obj.d), obj.d = -obj.RbarChol'\(obj.RbarChol\obj.d); end
                end
            end
        end % end of method: factor
        
        function str = double(obj)
            str = struct('class', 'MPCGPAD', 'sys', double(obj.sys), ...
                'N', obj.N, 'Q', obj.Q, 'R', obj.R, 'S', obj.S,...
                'q', obj.q, 'r', obj.r, 'QN', obj.QN, 'qN', obj.qN, ...
                'FN', obj.FN, 'gN', obj.gN, 'K', obj.K, 'D', obj.D, ...
                'M', obj.M, 'd', obj.d, 'L', obj.L, 'C', obj.C, ...
                's', obj.s, 'y', obj.y, 'nf', obj.nf, ...
                'is_q_zero', obj.is_q_zero, ...
                'is_r_zero', obj.is_r_zero, ...
                'is_S_zero', obj.is_S_zero);
        end % end of method: double
        
    end % end of all public methods
    
end %end of class: MPCGPAD

