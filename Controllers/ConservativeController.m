classdef ConservativeController < Controller
    % ConservativeController is a pseudo-controller to counteract
    % numerical losses
    %   This controller is only designed for (unactuated) conservative
    %   mechanical systems. ConservativeController counteracts energy
    %   losses in numerical integration schemes. Thus, enables the search
    %   for conservative limit cycles with fast but poor integrators. On
    %   the other hand, the root-finding-problem for a limit cycle becomes 
    %   quadratic with the additional decision variable xi.
    %
    %   Properties:
    %       nXi          number of control parameters
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
    %  see also: Controller
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 13-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        nXi       = 1       % number of free control parameters
        timeBased = false
    end
    methods (Access = public)
        %% - constructor -
        function obj = ConservativeController(model, lambda_E0)
            % constructor
            obj.instantData.M = zeros(model.nQ);
            obj.instantData.B = zeros(model.nQ,1);
            obj.xi            = lambda_E0;
        end
        %% - controller function, generates the input for the system -
        function [u,u_t,u_x,u_xi,u_T] = inputU(obj, ~, x, z, model, varargin)
            % u -> xDot = f(x) + u
            epsilon = model.getEpsilon();
            E_x     = model.Grad1.E_x(epsilon, x, z);
            E_xx    = model.Grad1.E_xx(epsilon, x, z);
            xi      = obj.xi;
            
            u       = E_x'*xi;
            u_x     = E_xx*xi;
            u_xi    = E_x';

            u_t     = [];
            u_T     = [];
        end
        
        function [tau, tau_t, tau_x, tau_xi, tau_T] = inputTau(obj, ~, ~, ~, model, varargin)
            % tau -> M*ddq = B*tau + ...
            nTau     = model.nTau;
            tau      = zeros(nTau,1);
            tau_x    = zeros(nTau,2*model.nQ);
            tau_xi   = zeros(nTau,obj.nXi);

            tau_t = [];
            tau_T = [];
        end
 
        
        function param = getFreeParameters(obj)
            % set xi
            param = obj.xi;
        end

        function setFreeParameters(obj, param)
            % get xi
            obj.xi = param;
        end

    end
end