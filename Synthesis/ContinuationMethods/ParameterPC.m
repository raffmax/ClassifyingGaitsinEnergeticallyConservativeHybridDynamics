classdef ParameterPC < PredictorCorrector
    % ParameterPC is a blueprint for a continuation method.
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
        function obj = ParameterPC(RFP,LC) 
            obj.RFP = RFP;
            obj.LC  = LC;
        end
    end
    methods
        function predictor(obj)
        end
        function corrector(obj)
        end
    end
end