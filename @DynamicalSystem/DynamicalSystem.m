classdef DynamicalSystem   
    properties       
        % Name of the system
        name = 'Dynamical System';
        % Note
        note = '';
    end % end of public properties
    
    
    properties(Access=protected)
         % Dimension of the state vector
        nx=0;
        % Number of inputs
        nu=0;
        % Number of inequality constraints
        nc=0;
        x=[];
    end % end of protected propertied
    
    methods
        
        function obj = setState(obj, x)
            obj.x=x;
        end
        
        function x = getState(obj)
            x = obj.x;
        end
        
        function [nx, nu, nc] = getDimensions(obj)
            nx = obj.nx;
            nu = obj.nu;
            nc = obj.nc;
        end
    end
end % end of class : DynamicalSystem