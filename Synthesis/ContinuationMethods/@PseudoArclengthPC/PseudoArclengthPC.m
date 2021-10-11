classdef PseudoArclengthPC < PredictorCorrector
    % PseudoArclengthPC is a predictor-corrector continuation method.
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
    % Last revision: 20-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       kappaTilde    % nominal contraction rate of corrector process
       deltaTilde    % nominal first corrector steplength
       alphaTilde    % nominal angel between two consecutive steps
       fixedStepSize % fixed or variable step size in predictor
       MaxIterPred   % maximum number of iterations in predictor
    end
    methods 
        %  - constructor - 
        function obj = PseudoArclengthPC(RFP,LC,varargin) 
            obj.constructorPC(RFP,LC,varargin{:});
            % add further PC settings
            % set up input parser
            default_kappa         = 1e-3;
            default_delta         = 1e-2;
            default_alpha         = pi/16;
            default_fixedStepSize = false;
            default_MaxIterPred   = 5;
            
            p  = inputParser;
            p.KeepUnmatched = true;
            
            addParameter(p, 'kappaTilde',    default_kappa         ,@isnumeric);
            addParameter(p, 'deltaTilde',    default_delta         ,@isnumeric);
            addParameter(p, 'alphaTilde',    default_alpha         ,@isnumeric);
            addParameter(p, 'fixedStepSize', default_fixedStepSize ,@islogical);
            addParameter(p, 'MaxIterPred',   default_MaxIterPred   ,@isnumeric);
            
            parse(p, varargin{:}); 
            
            obj.kappaTilde    = p.Results.kappaTilde;
            obj.deltaTilde    = p.Results.deltaTilde;
            obj.alphaTilde    = p.Results.alphaTilde;
            obj.fixedStepSize = p.Results.fixedStepSize;
            obj.MaxIterPred   = p.Results.MaxIterPred;
            
            
        end
    end

    
    methods       
        testfunctionsRoot(obj,LCBefore, LCAfter) % test for roots of various functionals
        
        function var = getConVar(obj,decVar,conPar)
            if obj.idxConPar > length(decVar)
                var = [decVar;conPar]; % conPar needs to be appended
            else
                var = decVar; % conPar already included
            end
        end
        
        function [H,HPrime_u] = getH(obj,varargin)
           % compute H(u) and its jacobian HPrime_u
           % copy solution struct
           u       = varargin{1};
           x       = u([1:(obj.idxConPar-1),obj.idxConPar+1:end]);
           lambda  = u(obj.idxConPar);
           calcJac = true; 
           if nargin == 3
               calcJac = varargin{2};
           end
           sol = obj.getSolStruct(x,lambda);
           
           switch obj.Type
               case 'E0'
                   if calcJac
                       [H,Jacobian] = getFandJacobian(obj.RFP,sol,x);
                       col          = zeros(size(x));
                       col(obj.RFP.idxCon.Additional(1)) = -1;
                       HPrime_u     = [Jacobian,col];
                   else
                       H        = getF(obj.RFP,sol,x);
                       HPrime_u = [];
                   end
               case 'epsilon'
                   obj.RFP.fixEpsilon(false);
                   if calcJac
                       [H,HPrime_u] = getFandJacobian(obj.RFP,sol,u);
                   else
                       H        = getF(obj.RFP,sol,u);
                       HPrime_u = [];
                   end
               case {'xDot_avg','E_avg'} 
                   if calcJac
                       [H,Jacobian] = getFandJacobian(obj.RFP,sol,x);
                       col          = zeros(size(x));
                       col(obj.RFP.idxCon.Additional(1)) = -sum(sol.tDomain);
                       HPrime_u     = [Jacobian,col];
                   else
                       H        = getF(obj.RFP,sol,x);
                       HPrime_u = [];
                   end
               case 'xi'
                   
               case 'T'
                   warning('this needs work!')
%                    indexT = length(sol.x0)+1;
                    u = [x;lambda];
                    [H,HPrime_u] = getFandJacobian(obj.RFP,sol,u);
%                    gradT = Jacobian(:,indexT);
%                    Jacobian(:,indexT) = [];
%                    HPrime_u     = [Jacobian, gradT];
               
               otherwise
           end

           function F = getF(objRFP,sol,X)
              fcn = @(X) objRFP.rootFunctionTDsingle(sol,X,objRFP.optimalControl);
              grad1 = objRFP.options.Grad1;
              
              objRFP.options.Grad1 = false;
              F = fcn(X);
              objRFP.options.Grad1 = grad1;
           end
           function [F,Jacobian] = getFandJacobian(objRFP,sol,X)
              fcn = @(X) objRFP.rootFunctionTDsingle(sol,X,objRFP.optimalControl);
              if objRFP.options.Grad1 
                  [F,Jacobian] = fcn(X);
              else
                  stepSize = objRFP.options.FiniteDifferenceStepSize;
                  [F,Jacobian] = getJacobianFD(fcn,X,stepSize);
              end
           end
            
        end
           
    end
end