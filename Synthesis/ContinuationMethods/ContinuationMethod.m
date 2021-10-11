classdef (Abstract) ContinuationMethod < matlab.mixin.Copyable
    % ContinuationMethod is a blueprint for a continuation method.
    %   This abstract model provides the important attributes and methods
    %   for a mechanical model homotopy.
    %
    %   Properties:  
    %
    %
    %   Methods:
    %  
    %
    % see also: RootFindingProblem, LimitCycle
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 15-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public) % --- References ---
        RFP  % root finding problem 
        LC   % limit cycle object of a limit cycle tree
    end
    
    properties (Access = public) 
        Tree
        Bifurcations
    end
    
%     properties (Access = public)
%         x      % free variables of problem (decision Variable)
%         lambda % continuation/homotopy parameter
%     end

    methods 
    end
end

