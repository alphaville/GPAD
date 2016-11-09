classdef LTISystem < DynamicalSystem
    %LTISystem is a dynamical system of the form: x+ = Ax + Bu + f which is
    %accompanied by a set of set and input constraints of the form:
    % - Box state constraints: xmin <= x <= xmax
    % - General state/input constraints: cmin <= Fx + Gu <= cmax
    % - Box input constraints: umin <= u <= umax
    %
    properties(Access=private)
        % System matrices: x+ = Ax + Bu + f
        A=[];
        B=[];
        f=[];
        % General inequality constraints: cmin <= Fx + Gu <= cmax
        F=[];
        G=[];
        cmin=[];
        cmax=[];
        % Input bounds: umin <= u <= umax
        umin=[];
        umax=[];
        % State bounds: xmin <= x <= xmax
        xmin=[];
        xmax=[];
        % Whether f=0:
        is_f_zero=false;
    end % End of public properties        
    
    methods
        
        function obj = LTISystem(varargin)
            %LTISystem can be constructed with one of the following ways:
            % - s = LTISystem();
            %    This will create an empty LTISystem object.
            % - s = LTISystem(other);
            %    This will create a copy of an existing LTISystem object.
            % - s = LTISystem(str); where str is a structure with the
            %    fields: A, B, f, F, G, cmin, cmax, xmin, xmax, umin, umax,
            %    name.
            % - s = LTISystem(A, B, f, F, G, cmin, cmax, xmin, ...
            %               xmax, umin, umax, name);
            %
            obj.note='LTI Dynamical System';
            if nargin==0
                %No arguments => Empty LTISystem object
                return;
            end
            if (nargin==1 && isa(varargin{1},'LTISystem'))
                args=varargin{1};
                %Copy another LTI System
                obj.name = args.name;
                obj.A = args.A;
                [obj.nx, obj.nu]=size(obj.B);
                obj.B = args.B;
                obj.f = args.f;
                obj.F = args.F;
                obj.G = args.G;
                obj.cmin = args.cmin;
                obj.cmax = args.cmax;
                obj.xmin = args.xmin;
                obj.xmax = args.xmax;
                obj.umin = args.umin;
                obj.umax = args.umax;
                obj.name = args.name;
                return;
            end
            if nargin>=2,
                obj.A = varargin{1};
                obj.B=varargin{2};
                obj.is_f_zero=true;
                [obj.nx, obj.nu]=size(obj.B);
                if (nargin==2), return; end;
                obj.f=varargin{3};
                if (isempty(obj.f) || norm(obj.f)<eps),
                    obj.is_f_zero=true;
                end
                if (nargin==3),
                    return;
                elseif (nargin<7)
                    error('Wrong number of input variables!');
                end;
                obj.F=varargin{4};
                obj.nc = size(obj.F, 1);
                obj.G=varargin{5};
                obj.cmin=varargin{6};
                obj.cmax=varargin{7};
                if (nargin==7), return; end;
                obj.xmin=varargin{8};
                if (nargin==8), return; end;
                obj.xmax=varargin{9};
                if (nargin==9), return; end;
                obj.umin=varargin{10};
                if (nargin==10), return; end;
                obj.umax=varargin{11};
                if (nargin==11), return; end;
                obj.name=varargin{12};
                if (nargin==12),
                    return;
                else
                    error('Too many input variables!');
                end;
            end
        end % end of constructor
        
        function display(obj)
            disp('LTI System : x+ = Ax + Bu');
            disp(['  ' num2str(obj.nx) ' states, ' num2str(obj.nu) ...
                ' inputs, ' num2str(obj.nc) ' constraints.']);
        end
        
        function obj = setDynamics(obj, A, B, f)
            obj.A = A;
            obj.B = B;
            if (varargin==4)
                obj.f = f;
            end
        end
        
        function [A, B, f] = getDynamics(obj)
            A = obj.A;
            B = obj.B;
            f = obj.f;
        end
        
        function obj = stateupdate(obj, u)
            %STATEUPDATE computes the next state for a given input action u.
            obj.x = obj.A * obj.x + obj.B * u;
            if ~obj.is_f_zero, obj.x = obj.x + obj.f; end
        end
        
        function obj = setConstraints(obj, F, G, cmin, cmax, umin, umax, ...
                xmin, xmax)
            %SETCONSTRAINTS is a setter method for the constraints of the
            %LTI System:
            %
            % (x,u): cmin <= Fx + Gu <= cmax
            %        xmin <= x <= xmax
            %        umin <= u <= umax
            obj.F=F;
            obj.G=G;
            obj.cmin=cmin;
            obj.cmax=cmax;
            obj.umin=umin;
            obj.umax=umax;
            obj.xmin=xmin;
            obj.xmax=xmax;
            obj.nc = size(F,1);
        end
        
        function [F, G, cmin, cmax, umin, umax, xmin, xmax] = ...
                getConstraints(obj)
            %GETCONSTRAINTS returns the matrices of the contraints of the
            %LTI system.
            F = obj.F;
            G=obj.G;
            cmin=obj.cmin;
            cmax=obj.cmax;
            umin=obj.umin;
            umax=obj.umax;
            xmin=obj.xmin;
            xmax=obj.xmax;
        end
        
        function ifz = isfzero(obj)
            ifz = obj.is_f_zero;
        end
        
        function str = double(obj)
            str = struct('class', 'LTI', 'A', obj.A, 'B', obj.B, ...
                'f', obj.f, 'F', obj.F, 'G', obj.G, 'cmin', obj.cmin, ...
                'cmax', obj.cmax, 'umin', obj.umin, 'umax', obj.umax, ...
                'xmin', obj.xmin, 'xmax', obj.xmax, 'is_f_zero', obj.is_f_zero, ...
                'nx', obj.nx, 'nu', obj.nu, 'nc', obj.nc, 'name', obj.name, ...
                'x', obj.x);
        end
        
        %TODO: Method #createRandom(nx,nu,nc)
    end % end of public methods
    
end