classdef (Abstract) Controller < handle
    % Controller is a blueprint for different control structures
    %   This abstract model provides the important attributes and methods
    %   for a model homotopy controller.
    %
    %   Properties:
    %       nXi          number of free control parameters
    %       model        reference to mechanical model
    %       xi           all control parameters
    %       instantData  instantaneous data ( M(q(t)), B(t), ...) at time t     
    %
    %   Methods:
    %       inputU       computes general input u in xDot = f(x) + u
    %       inputTau     computes torque input tau in M*ddq = B*tau + ...
    %       update       updates M and B at time instance t
    %
    %       getFreeParameters   get free ctrl parameter (decision variable)
    %       setFreeParameters   set free ctrl parameter (decision variable)
    %      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 13-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Abstract)
        nXi           % number of free control parameters
        timeBased     % specify if controller is time dependent
    end
    
    properties (Access = protected)
        model         % reference to mechanical model
        xi            % all control parameters
        instantData   % instantaneous data ( M(q(t)), B(t), ...) at time t
    end
    
    
    methods (Abstract) %% - controller function, generates the input for the system -
        % u -> xDot = f(x) + u
        [u, u_t, u_x, u_xi, u_T]         = inputU(obj, t, x, z, model, varargin)
        % tau -> M*ddq = B*tau + ...
        [tau, u_t, tau_x, tau_xi, tau_T] = inputTau(obj, t, x, z, model, varargin)
       
        param = getFreeParameters(obj) % get free ctrl parameter (decision variable)
        setFreeParameters(obj, param)  % set free ctrl parameter (decision variable)
        
    end
    
    methods %% - update -
        function update(obj, B, M)
            % update overwrites instantaneous B and M
            obj.instantData.B = B;
            obj.instantData.M = M;
        end
        function B = getB(obj)
            B = obj.instantData.B;
        end
        function M = getM(obj)
            M = obj.instantData.M;
        end
    end
end