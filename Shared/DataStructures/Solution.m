classdef Solution
    % Solution represents a container class for results.
    %   Solution stores results of any kind for the purpose of finding
    %   limit cycles. This is not a handle class!
    %   The following properties are required to uniquelly
    %   obtain a system's behaviour by simulation in RootFindindingProblem
    %   or LimitCycle.
    %
    %   Properties:  
    %      x0            number of (physical) constraints
    %      tDomain       number of generalized coordinates
    %      xi            controll parameters
    %      epsilon       homotopy parameters
    %
    %      rfpData       data struct of root-finding-problem, with fields, e.g decVar
    %      param         solution parameters, e.g. E0, xDot_avg, E_avg
    %      dynamics      data struct with fields MonodromyMatrix and f0
    %
    %
    %   Methods:
    %
    %   Example:
    %       x0      = [1;0];
    %       tDomain = [2,1,2];
    %       xi      = 0;
    %       epsilon = 0;
    %       sol = Solution(x0,tDomain,xi,epsilon);
    %
    %   see also RootFindindingProblem, LimitCycle
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 05-May-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
        x0            % number of (physical) constraints
        tDomain       % number of generalized coordinates
        xi            % controll parameters
        epsilon       % homotopy parameters
    end
    
    properties (Access = public)
        rfpData       % data struct of root-finding-problem, with fields, e.g decVar
        param         % solution parameters, e.g. E0, xDot_avg, E_avg
	    dynamics      % data struct with fields MonodromyMatrix and f0
    end
    
    methods
        function sol = Solution(x0,tDomain,xi,epsilon,varargin)
            %Constructor

            p = inputParser;
            
            addRequired(p,'x0'     ,@isnumeric)
            addRequired(p,'tDomain',@isnumeric)
            addRequired(p,'xi'     ,@isnumeric)
            addRequired(p,'epsilon',@isnumeric)

            parse(p,x0,tDomain,xi,epsilon);
            
            for k=1:numel(p.Parameters)
                field = p.Parameters{k};
                sol.(field) = p.Results.(field);
            end
            
            if nargin>4
                oldSolution = varargin{:};
                if ~isequal(class(oldSolution),'Solution')
                    error('Extra argument must be a solution object!')
                end
                % copy remaining properties
                props = properties(oldSolution);
                for iprop = 1:length(props)
                    thisprop = props{iprop};
                    if isempty(sol.(thisprop))
                        sol.(thisprop) = oldSolution.(thisprop);
                    end
                end
            else
                sol.rfpData.multiplier = [];
            end
               
        end
    end
end