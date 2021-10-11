classdef (Abstract) PredictorCorrector < ContinuationMethod
    % PredictorCorrector is a blueprint for a continuation method.
    %   This abstract model provides the important attributes and methods
    %   for a mechanical model homotopy.
    %
    %   Properties:  
    %
    %
    %   Methods:
    %  
    %
    % see also: ContinuationMethod, RootFindingProblem
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 15-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        StepSize      % nominal step size (h) in Corrector
        StepTolerance % checks if step size is too small -> stop continuation
        MaxIterations % max. number of PC steps
        Direction     % direction +1/-1 (e.g. in tangent space orientation) 
        Type          % specifies type of continuation parameter
        dispIter      % whether or not to display iterations
        dispBifurcation % display found bifurcations; see in check4bifurcation.m
        idxConPar     % index of continuation Parameter
    end
    
    methods (Access = protected)
        function constructorPC(obj,RFP,LC,varargin)
            % general constructor for PC methods
            obj.RFP = RFP;
            obj.LC  = LC;
            % set up input parser
            default_StepSize      = 0.05; 
            default_MaxIterations = 1000;
            default_StepTolerance = 1e-5;
            default_Direction     = +1;
            default_dispIter      = false;
            default_dispBif       = true;
            default_Tree          = [1,1];
            
            p = inputParser;
            p.KeepUnmatched =  true; % allows to set parameters with  a parser e.g. in the PseudoArclengthPC constructor
            addParameter(p, 'StepSize',       default_StepSize       ,@isnumeric);
            addParameter(p, 'MaxIterations',  default_MaxIterations  ,@isnumeric);
            addParameter(p, 'StepTolerance',  default_StepTolerance  ,@isnumeric);
            addParameter(p, 'Direction',      default_Direction      ,@isnumeric);
            addParameter(p, 'dispIter',       default_dispIter       ,@islogical);
            addParameter(p, 'dispBifurcation',default_dispBif        ,@islogical);
            addParameter(p, 'Tree',           default_Tree           ,@isnumeric);
            
            
            parse(p, varargin{:}); 
            
            obj.StepSize      = p.Results.StepSize;
            obj.MaxIterations = p.Results.MaxIterations;
            obj.StepTolerance = p.Results.StepTolerance;
            obj.Direction     = p.Results.Direction;
            obj.dispIter      = p.Results.dispIter;
            obj.dispBifurcation=p.Results.dispBifurcation;
            obj.Tree          =p.Results.Tree;
            
            if isfield(obj.LC.Sol.rfpData, 'direction') 
                if ~isempty(obj.LC.Sol.rfpData.direction)
                    obj.Direction = obj.LC.Sol.rfpData.direction;
                end
            end
        end
    end
    
    methods
        function [exitflag,stepsPC] = runContinuation(obj,type,varargin)
           obj.Type = type;
           % get value of continuation parameter (conPar) from current Limit-Cycle
           obj.idxConPar     = obj.getIdxConPar();
           switch type
               case 'E0'
                   conPar = obj.LC.model.SystemEnergy(obj.LC.Sol.x0,obj.LC.sequence.FlowMapOrder(:,1));
               case 'epsilon'
                   conPar = obj.LC.Sol.epsilon;
               case 'xDot_avg'
                   conPar = obj.LC.Sol.param.xDot_avg;
               case 'xi'
                   
               case 'E_avg'
                   conPar = obj.LC.Sol.param.E_avg;   
               case 'T'
                   conPar = obj.LC.Sol.tDomain;           
               otherwise
           end
           disp(['Run PC-Continuation in ',type,' ...'])
           if nargin >2
               target = varargin{:};
               [exitflag,stepsPC] = runPC(obj,conPar,target);
           else
               [exitflag,stepsPC] = runPC(obj,conPar);
           end
           if exitflag == 1
               disp('... terminated successfully with exitflag 1.')
           else
               disp(['... Continuation failed! See exitflag ',num2str(exitflag),'.'])
           end
        end
        
        function sol = getSolStruct(obj,varargin)
           if nargin < 3
               var    = varargin{1};
               decVar = var(1:end-1);
               conPar = var(end);
           else
               decVar = varargin{1};
               conPar = varargin{2};
           end
            
            switch obj.Type
               case 'E0'
                   sol                = obj.LC.Sol;
                   sol.rfpData.decVar = decVar;
                   sol.param.E0       = conPar; 
               case 'epsilon'
                   % update epsilon in solLC
                   solLC              = obj.LC.Sol;
                   sol                = Solution(solLC.x0,solLC.tDomain,solLC.xi,conPar,solLC);
                   sol.rfpData.decVar = decVar;
                   obj.LC.model.setEpsilon(conPar);
               case 'xDot_avg'
                   sol                = obj.LC.Sol;
                   sol.rfpData.decVar = decVar;
                   sol.param.xDot_avg = conPar; 
                   
               case 'xi'
                   
               case 'E_avg'
                   sol                = obj.LC.Sol;
                   sol.rfpData.decVar = decVar;
                   sol.param.E_avg    = conPar; 
                   
                case 'T'
                   sol                = obj.LC.Sol;
                   sol.rfpData.decVar = decVar;
                   sol.param.T        = sol.tDomain;
               otherwise
           end
        end
        
        function [conPar,name] = getConPar(obj,sol)
            switch obj.Type
               case 'E0'
                   conPar = sol.param.E0;
                   name = 'E0';
               case 'epsilon'
                   conPar = sol.epsilon; 
                   name = 'epsilon';
               case 'xDot_avg'
                   conPar = sol.param.xDot_avg;
                   name = 'xDot_avg';
               case 'xi'

               case 'E_avg'
                   conPar = sol.param.E_avg;
                   name = 'E_avg';
                case 'T'
                   conPar = sol.param.T;
                   name = 'T';
               otherwise
           end
        end
        
        function dist = getSolDistance(obj, solBase, solTarget)
            % compute (directional) steplenght between solution solBase
            % and solTarget. the sign of dist corresponds to the search
            % direction needed to get from solBase to solTarget
            
            conParBase   = obj.getConPar(solBase);
            conParTarget = obj.getConPar(solTarget);
            
            uBase    = [solBase.rfpData.decVar; conParBase];
            uTarget  = [solTarget.rfpData.decVar; conParTarget];
            
            dist     = norm(uTarget-uBase) * solBase.rfpData.direction; 
        end
        
        function newLC = addLC(obj,sol,d,varargin)
            % input: sol   - solution for which a LC is added
            %        d     - wether to add it after/next (+1) the current
            %                obj.LC or before/prev (-1)
            % optional input: 
            %        {1} a LC to which the newLC is added
            
            switch obj.Type
               case 'E0'
                   E0     = sol.param.E0;
                   name   = ['E0=',num2str(E0)];
               case 'epsilon'
                   epsilon = sol.epsilon;
                   name    = ['eps=',num2str(epsilon)]; 
               case 'xDot_avg'
                   xDot_avg = sol.param.xDot_avg;
                   name     = ['xDot_avg=',num2str(xDot_avg)];
                   
               case 'xi'
                   
               case 'E_avg'
                   E_avg  = sol.param.E_avg;
                   name   = ['E_avg=',num2str(E_avg)]; 
                   
                case 'T'
                   T = sol.param.T;
                   name    = ['T=',num2str(T)]; 
               otherwise
            end    
            % create new limit cycle
            newLC  = LimitCycle(name,sol,obj.RFP);
            
            % next to which LC?
            if nargin < 4
                LC = obj.LC;
            else 
                LC = varargin{1};
            end
            
            % Prev or Next?
            if d == 1
                insertAfter(newLC, LC);
            else
                insertBefore(newLC, LC);
            end
        end
                
        
        function idxConPar = getIdxConPar(obj)
            switch obj.Type
               case 'E0'
                   idxConPar = obj.RFP.nDecVar + 1; %% needs to be checked!!
               case 'epsilon'
                   %idxConPar = obj.RFP.idx.Eps; 
                   idxConPar = obj.RFP.nDecVar+1;
               case 'x0'
                   
               case 'xDot_avg'
                  idxConPar = obj.RFP.nDecVar+1;

               case 'xi'

               case 'E_avg'
                   idxConPar = obj.RFP.nDecVar + 1; %% needs to be checked!!
               case 'T'
                   idxConPar = obj.RFP.idx.Time;
               otherwise
                   
           end
        end
        
        function [rootFcn, typeTestFcn, options] = getRootFindingFunctional(obj, type, varargin)
           % returns a functional for which a root is wanted; 
           % 
           % input: type (string) - 'detJacAug', 'eigJacAug'
           % output: rootFcn - function handle for root fcn with a LimitCycle as input 
           %         typeTestFcn - corresponds to the value assigned in
           %                       testfunctionsRoot if a root is detected
           %         options - optional output e.g. an alternative stopping
           %                   criteria 
           
           options = struct(); 
           
            switch type
                case 'detJacAug'
                    rootFcn = @(LC) det(LC.Sol.rfpData.jacobianAug)*LC.Sol.rfpData.direction;
                    options.stopCrit = @(LC) eigs(LC.Sol.rfpData.jacobianAug, 1, 'smallestabs'); % smallest absolute eig value
                    typeTestFcn = 1;
                case 'tDomain'
                    if ~isempty(varargin)
                        error('The timedomain for which a root should be found must be passed as varargin!')
                    else
                        tIdx = varargin{1};
                    end
                    rootFcn =  @(LC) LC.Sol.tDomain(tIdx);
                otherwise
                    error("Type of functional for which a root is wanted is not implemented (yet).");
            end
            
        end
        
    end
    
    methods (Abstract) % PC-Continuation in continuation parameter
        [exitflag,stepsPC] = runPC(obj,conPar,target)
    end
    
    
    methods (Abstract)
        predictor(obj)
        corrector(obj)
    end
end
