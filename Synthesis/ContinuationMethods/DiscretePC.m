classdef DiscretePC < PredictorCorrector
    % DiscretePC is a predictor-corrector continuation method.
    %   This abstract model provides the important attributes and methods
    %   for a mechanical model homotopy.
    %
    %   Properties:  
    %
    %
    %   Methods:
    %  
    %
    % see also: ContinuationMethod, RootFindingProblem, LimitCycle
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 15-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods 
        %  - constructor - 
        function obj = DiscretePC(RFP,LC,varargin) 
            obj.constructorPC(RFP,LC,varargin{:});
        end
    end
    methods % PC-Continuation in continuation parameter
        function [exitflag,stepsPC] = runPC(obj,conPar,conPar_target)
            if nargin <3
                conPar_target = inf;
            end
            % epsilon is not a free variable in RFP
            obj.RFP.fixEpsilon(true);
            
            % get continuation direction
            d    = obj.Direction;
            iter = 0;
            
            stepsPC = d*obj.StepSize*ones(obj.MaxIterations,1);
           
            if obj.dispIter
                disp('Iteration:');
            end
            
            while d*conPar <= d*conPar_target && iter <= obj.MaxIterations
                if obj.dispIter
                    fprintf('%d ', iter);
                    if mod(iter, 50)==49 % print each 50 iterations in a new line
                        fprintf('\n');
                    end
                end 
                
                % get decison variable of current Limit Cycle
                LCsol  = obj.LC.Sol;
                decVar = obj.RFP.getDecisionVariable(LCsol);
                %% Predictor
                [decVarP,conParP] = obj.predictor(decVar,conPar);
                solP = obj.getSolStruct(decVarP,conParP);
                %% Corrector
                [solC,exitflagC] = obj.corrector(solP);
                
                %% add Limit Cycle
                if exitflagC > 0
                    % create and add new limit cycle
                    newLC  = obj.addLC(solC,d);
                    % update Limit Cycle and conPar
                    obj.LC = newLC;
                    conPar = obj.getConPar(obj.LC.Sol);
                else
                    warning('The corrector-step was not successful!')
                    break
                end
                
                % update
                iter = iter+1;
            end
            if iter < obj.MaxIterations
                exitflag = exitflagC;
            else
                warning('Continuation stopped prematurally')
                exitflag = 0;
            end
            % resize PC_steps
            stepsPC = stepsPC(1:iter-1); % -1 since iter[1]=0
            
            if obj.dispIter
                fprintf('\n') % end output of iterations
            end
            
        end
    end
    
    methods % implement Predictor & Corrector
        function [xP,lambdaP] = predictor(obj,x,lambda)
            xP      = x;
            lambdaP = lambda + obj.Direction*obj.StepSize;
        end
        function [solC,exitflagC] = corrector(obj,solP)
            [solC,~,exitflagC] = obj.RFP.solveNewtonRFP(solP);
        end
    end
end