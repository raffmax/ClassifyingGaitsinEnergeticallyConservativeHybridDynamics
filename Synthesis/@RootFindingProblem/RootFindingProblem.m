classdef RootFindingProblem < matlab.mixin.Copyable
    % RootFindingProblem defines the general structure of the root finding problem
    %   This 
    %
    %   Properties:  
    %
    %
    %   Methods:
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 30-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private) % --- Class-References ---
        model          % reference to mechanical model                                     
        controller     % reference to used controller
        simulation     % reference to simulation object                                                      
    end
    
    properties (SetAccess = private) % --- Model dependent props ---
        freeStates     % states that are part of decision variables
        fixedStates    % states that are not part of decision variables
        periodicStates % states on which we enforce periodicity with period T
        periodicStatesHomotopy % states that start being fixed (epsilon=0) and transition to periodicity (epsilon=1)

        sequence %Predefined Footfall pattern                                                         

        flipOperator % exploit symmetry or use coordinate transformation
        % - active implicit equality constraint (e.g. fix initial system energy)
        fixedParameter
        % - fix epsilon
        fixedEpsilon
        % - functionals for implicit equality constraints or costs
        functionals
        % - define optimization Problem in Control-Parameters
        optimalControl
        useLagrangeMultipliers
    end

    properties % --- Solver/Implementation dependent props ---
        integrationType % event-driven (ED), time-driven (TD), 
        shootingType    % single shooting (single), multiple shooting (multiple)
        nDecVar         % number of Decision Variables
        nCon            % number of Constraints
        nInt            % number of additional integrals
        idx             % struct indexing decVar for free-state, time and param 
        idxCon          % struct indexing equality constraints
        nState          % number of total states
        nFreeState      % number of free states
        nTime           % number of domain time
        nCtrl           % number of control parameters (not including epsilon)
        options         % tolerances, gradients, ...
    end
    


    methods 
        %  - constructor - 
        function obj = RootFindingProblem(model,sequence,varargin) 
            if nargin < 2
                error('To instantiate RootFinding() the minimum amount of arguments must be model and sequence!')
            end
            
            % check if model is of type MechanicalModel
            if ~strcmp(model.type,'MechanicalModel')
                disp(pogoStick)
                error('model must be of type MechanicalModel!')
            else
                obj.model    = model;
            end
            
            if size(sequence.FlowMapOrder,1) ~= model.nZ && model.nZ>0
                error(['dimension of discrete states in sequence must be ',num2str(model.nZ),'!'])
            else
                obj.sequence = sequence;
            end
                
            obj.generalSetup(varargin{:});
        end
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods (Access = private)
        generalSetup(obj,varargin) % setup in constructor
        nDecVar = getNumDecVar(obj)
        function updateNumDecVar(obj)
            obj.nDecVar = obj.getNumDecVar();
        end 
        nCon    = getNumConstraints(obj)
        nInt    = getNumIntegrands(obj)
        nTime   = getNumDomainTime(obj)
    end
    
    methods (Access = public) % get and set functions
        Lxx = getHessian(obj,rootFun,x,lambda)
        setOptimalControl(obj,optLogical)
        setFunctional(obj,idx,type,varargin)
        removeFunctional(obj,idx)
        setSolverOptions(obj,options)
        fixEpsilon(obj,epsLogical)
        fixParameter(obj,paramLogical)
        decVar  = getDecisionVariable(obj,x0,tDomain,epsilon)
        [sol,output] = getOutputSolvers(obj,init,decVar,output,jacobian,exitflag)
        function trajectory = getTrajectory(obj,x0,tDomain)
            trajectory = obj.simulation.getTrajectory(x0,tDomain);
        end
    end
    
    methods (Access = public) % root functions
        % compute rootFunction by time-driven single shooting
        [fun,jacobian,objective] = rootFunctionTDsingle(obj,init,decVar,computeLxx)
    end

    methods (Access = public) % - Root-Finding-Methods
        % executable solvers and functions
        [sol,fval,exitflag,output,jacobian] = fsolveRFP(obj,init)
        [sol,fval,exitflag,output,jacobian] = solveNewtonRFP(obj,init)
    end

    methods (Access = private)
        % additional methods for executable solvers
        options   = setOptionsfsolve(obj)
        [xE,intE] = getGradientsFromSimData(obj,simData)
    end

end

