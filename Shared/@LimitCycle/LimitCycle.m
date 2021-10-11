classdef LimitCycle < matlab.mixin.Copyable
    % LimitCycle Represents a periodic solution.
    %   A limit cycle is a periodic solution of MechanicalModel and can be
    %   related to other limit cycles to form a connected tree.
    %
    %   Properties:
    %
    %   RigidBody Methods:
    %
    %   See also RootFindingProblem, MechanicalModel, Controller,
    %   Simulation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 14-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private) % --- References ---
        model      % reference to MechanicalModel
        controller % reference to Controller
        simulation % reference to Simulation object
    end
    
    properties (SetAccess = private)
        Name                              % name of limit cycle
        Label
        
        Bifurcation = []                  % is this LC a root (e.g. bifurcation)
        BeforeBif   = []                  % is this LC located before a root
        AfterBif    = []                  % is this LC located after a root
        
        Prev        = LimitCycle.empty    % handle of previous limit cycle
        Next        = LimitCycle.empty    % handle of next limit cycle
        Sol                               % reference to Solution object
        
        sequence                          % predefined footfall pattern / contact sequence
    end
    
    
    properties (SetAccess = private) % required properties for simulation
        flipOperator    % here identity!
        nState          % number of total states
        nTime           % number of domain time
        nCtrl           % number of control parameters (not including epsilon)
        nInt            % number of additional integral-functionals
    end
    
    methods
        %  - constructor -
        function obj = LimitCycle(name,sol,varargin)
            %LimitCycle Constructor
            %  varargin can be:
            %           - RootFindingProblem object
            %           - model,controller,simulation,sequence
            if nargin <5
                rfp = varargin{1}; % possibly object of RootFindingProblem
                if ~isa(rfp, 'RootFindingProblem')
                    error('varargin must be of type RootFindingProblem')
                end
                
                if isprop(rfp,'flipOperator')
                    obj.flipOperator = rfp.flipOperator;
                else
                    obj.flipOperator = eye(obj.nState);
                end
            else
                p = inputParser;
                % cfg
                addRequired(p, 'model');
                addRequired(p, 'controller');
                addRequired(p, 'simulation');
                addRequired(p, 'sequence');
                
                parse(p, varargin{:});
                rfp = p.Results;
                
                obj.flipOperator = rfp.flipOperator;
            end
            obj.model      = rfp.model;
            obj.controller = rfp.controller;
            % set new simulation tolerance
            obj.simulation = Simulation('odeOptions',odeset('RelTol',1e-5,'AbsTol',1e-7));
            obj.sequence   = rfp.sequence; %TO-DO: careful if sequence exploits symmetry
            obj.simulation.setSystemRef(obj); % set new reference in simulation object
            
            obj.Name        = name;
            obj.Label       = [1,1];
            obj.Sol         = sol;
            
            % set Simulation properties
            obj.nState       = 2*obj.model.nQ;
            obj.nCtrl        = length(sol.xi);
            obj.nInt         = 0;
            obj.nTime        = length(sol.tDomain);
        end
        
        conData = getContinuationData(LC)
        trajectories = getContinuationTrajectories(LC,type,conParRange,conParSize,tSize,direction)
        
        function updateLCSol(LC,solIN)
            LC.Sol.rfpData  = solIN.rfpData;
            LC.Sol.param    = solIN.param;
            LC.Sol.dynamics = solIN.dynamics;
        end
        
        function setController(LC,controller)
            LC.controller = controller;
        end
        
        function insertAfter(newLC, prevLC)
            removeLC(newLC);
            newLC.Next = prevLC.Next;
            newLC.Prev = prevLC;
            if ~isempty(prevLC.Next)
                prevLC.Next.Prev = newLC;
            end
            prevLC.Next = newLC;
        end
        
        function insertBefore(newLC, nextLC)
            removeLC(newLC);
            newLC.Next = nextLC;
            newLC.Prev = nextLC.Prev;
            if ~isempty(nextLC.Prev)
                nextLC.Prev.Next = newLC;
            end
            nextLC.Prev = newLC;
        end
        
        function removeLC(LimitCycle)
            if ~isscalar(LimitCycle)
                error('Limit Cycles must be scalar')
            end
            prevLC = LimitCycle.Prev;
            nextLC = LimitCycle.Next;
            if ~isempty(prevLC)
                prevLC.Next = nextLC;
            end
            if ~isempty(nextLC)
                nextLC.Prev = prevLC;
            end
            LimitCycle.Next = LimitCycle.empty;
            LimitCycle.Prev = LimitCycle.empty;
        end
        
        function disconnectLC(LimitCycle)
            if ~isscalar(LimitCycle)
                error('Limit Cycles must be scalar')
            end
            LimitCycle.Next = LimitCycle.empty;
            LimitCycle.Prev = LimitCycle.empty;
        end
        
        function clearList(LC)
            % Clear the list before
            % clearing list variable
            prev = LC.Prev;
            next = LC.Next;
            removeLC(LC)
            while ~isempty(next)
                LC = next;
                next = LC.Next;
                removeLC(LC);
            end
            while ~isempty(prev)
                LC = prev;
                prev = LC.Prev;
                removeLC(LC)
            end
        end
        
        function [LC, flag] = getNextLC(obj)
            if length(obj.Next) == 1 && ~isnan(obj.Next.Label(1))
                LC      = obj.Next;
                flag    = 1;
            elseif isempty(obj.Next)
                LC = obj;
                flag = 0;
            else
                disp('Next limit cycle is a bifurcation point. Choose next limit cycle manually.')
                LC = obj;
                flag = 0;
            end
        end
         
        function [LC, flag] = getPrevLC(obj)
            if length(obj.Prev) == 1 && ~isnan(obj.Prev.Label(1))
                LC      = obj.Prev;
                flag    = 1;
            elseif isempty(obj.Prev)
                LC = obj;
                flag = 0;
            else
                disp('Previous limit cycle is a bifurcation point. Choose previous limit cycle manually.')
                LC = obj;
                flag = 0;        
            end
        end
        
        function LC = getLCstart(obj)
            if ~isempty(obj.Prev)
                [LC,flag] = obj.getPrevLC();
                if flag == 0
                    LC = obj;
                    return
                end
                LC = LC.getLCstart();
            else
                LC = obj;
            end
        end
        
        function LC = getLCend(obj)
            if ~isempty(obj.Next) && ~isnan(obj.Next.Label(1))
                [LC,flag] = obj.getNextLC();
                if flag == 0
                    LC = obj;
                end
                LC = LC.getLCend();
            else
                LC = obj;
            end
        end
        
        function depth = treeDepth(obj)
            depth = 0;
            LC = obj;
            flag = 1;
            while ~isempty(LC.Next) && flag ~= 0
                depth = depth + 1;
                [LC,flag] = LC.getNextLC();
            end
        end
        
        function addPrev(obj,varargin)
            for i = 1:length(varargin)
                obj.Prev(end+1) = varargin{i};
                varargin{i}.Next(end+1) = obj;
            end
        end
        
        function addNext(obj,varargin)
            for i = 1:length(varargin)
                obj.Next(end+1) = varargin{i};
                varargin{i}.Prev(end+1) = obj;
            end
        end
        
        function varargout = findNextBifurcation(obj)
            % input:  obj - LC to start search in the direction of LC.Next
            % output: if bif is already approximated: [LCBif]
            %         else [LCBefore; LCAfter]
            
            LC = obj; % obj = LCstart
            flag = 1;
            while isempty(LC.BeforeBif) && isempty(LC.AfterBif) &&  ~isempty(LC.Next) && flag ~= 0
                [LC,flag] = LC.getNextLC();
            end
            if ~isempty(LC.Next) && isempty(LC.BeforeBif) && flag ~= 0% LC is LC.Next
                varargout{1} = LC.Next;
                varargout{2} = LC;
            elseif ~isempty(LC.Next) && isempty(LC.AfterBif) && flag ~= 0% LC is LC.Prev
                varargout{1} = LC;
                varargout{2} = LC.Next;
            elseif ~isempty(LC.Bifurcation)
                varargout{1} = LC;
                varargout{2} = LC.empty();
                disp('Bifurcation Point reached. Switch to desired branch manually.')
            elseif ~isempty(LC.Next) && ~isempty(LC.Next.Bifurcation) % Bifurcation is already approximated
                varargout{1} = LC.Next;
                varargout{2} = LC.empty();
                disp('Bifurcation Point reached. Switch to desired branch manually.')
            end
            if isempty(LC.BeforeBif) && isempty(LC.AfterBif) && isempty(LC.Next)
                varargout{1} = LC.empty(); % no more bif in this LC tree
                varargout{2} = LC.empty();
            end
        end
        
        function trajectory = getTrajectory(obj,dt)
            obj.controller.setFreeParameters(obj.Sol.xi);
            obj.model.setEpsilon(obj.Sol.epsilon);
            if nargin < 2
                dt = [];
            end
            trajectory = obj.simulation.getTrajectory(obj.Sol.x0,obj.Sol.tDomain,dt);
        end
        
    end
    
    
    methods (Access = public) % set and clear functions for Bifurcations
        function setBifurcation(obj, type, varargin)
            % save as struct in obj.Bifurcation such that it is possible to
            % add optional details e.g. the dimension of the kernel
            if isempty(obj.Bifurcation)
                obj.Bifurcation = struct('type', type); % init struct
            else
                obj.Bifurcation(end+1).type = type; % add to struct
            end
            if ~isempty(varargin)
                obj.Bifurcation(end).details = varargin{1}; % add optional details 
            end
%             obj.Bifurcation = [obj.Bifurcation, type];
        end
        
        function setBeforeBifurcation(obj, type)
            obj.BeforeBif = [obj.BeforeBif, type];
        end
        
        function setAfterBifurcation(obj, type, varargin)
            obj.AfterBif = [obj.AfterBif, type];
        end
        
        function clearBeforeBifurcation(obj, varargin)
            if nargin>1
                type = varargin{1};
                obj.BeforeBif(obj.BeforeBif==type) = [];
            else
                obj.BeforeBif = [];
            end
        end
        function clearAfterBifurcation(obj, varargin)
            if nargin>1
                type = varargin{1};
                obj.AfterBif(obj.AfterBif==type) = [];
            else
                obj.AfterBif = [];
            end
        end
        
        function setLabel(obj, label)
           obj.Label = label; 
        end
        
        function Tree = updateLabel(obj, Tree, LC)
            if ~isempty(obj.AfterBif)
                idx = find(Tree(:,1)==obj.Label(1,1));
                newLabel = [obj.Label(1,1), max(Tree(idx,2))+1];
                obj.setLabel(newLabel);
                Tree = [Tree;newLabel];
            elseif ~isempty(obj.Bifurcation)
                obj.Label = NaN;
            elseif nargin > 2
                branchIdx = max(Tree(:,1))+1;
                idx = find(Tree(:,1)==branchIdx);
                if ~isempty(idx)
                    newLabel = [branchIdx, max(Tree(idx,2))+1];
                else
                    newLabel = [branchIdx, 1];
                end
                
                obj.setLabel(newLabel);
                LC.setLabel(newLabel+[0 1]);
                Tree = [Tree;newLabel;newLabel+[0 1]];
            end
        end
    end
    
    methods (Access = private)
        function delete(LC)
            clearList(LC)
        end
    end
end