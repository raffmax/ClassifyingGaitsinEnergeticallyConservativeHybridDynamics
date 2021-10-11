classdef Simulation < matlab.mixin.Copyable 
    properties (SetAccess = private) % --- References ---
        systemRef % reference to RootFindingProblem or LimitCycle
    end
    
    properties (Access = private)
        odeOptions
        odeSolver
    end
     
    methods % --- Constructor
        function obj = Simulation(varargin)
            % default values for the configuration
            default_odeOptions  = odeset('RelTol',1e-7,'AbsTol',1e-9);
            default_odeSolver   = @ode45;

            p = inputParser;
            % cfg 
            addParameter(p, 'odeOptions',    default_odeOptions);
            addParameter(p, 'odeSolver',     default_odeSolver);

            parse(p, varargin{:}); 

            obj.odeOptions        = p.Results.odeOptions;
            obj.odeSolver         = p.Results.odeSolver;
        end
        
    end
    
    methods (Access = public) % get and set methods
        function setSystemRef(obj,newSystemRef)
            obj.systemRef = newSystemRef;    
        end
        
        monodromy = computeMonodromyMatrix(obj,simData,eventGradients)
        
        function trajectory = getTrajectory(obj,x0,tDomain,dt,xi,epsilon,timeDriven)
            if nargin < 4
                dt = [];
            end
            if nargin > 4
                obj.systemRef.controller.setFreeParameters(xi);
                if nargin==5
                    obj.systemRef.model.setEpsilon(epsilon);
                end
                if nargin < 6
                    timeDriven = true;
                end
            else
                timeDriven = true;
            end
            T = sum(tDomain);
            if timeDriven
                % time driven
                [~,~,trajectory] = obj.singleShootingTD(x0,tDomain,false,false,true);
            else
                % event driven
                [~,~,trajectory] = obj.singleShootingED(x0,T,false,false,true);
            end

            % interpolate
            iDomain_del = []; % domains/rows that will be without information
            if ~isempty(dt)
                for iDomain = 1:length(trajectory)
                    tData = trajectory(iDomain).t;
                    t     = linspace(tData(1),tData(end),round((tData(end)-tData(1))/dt));
                    if ~isempty(t)
                        trajectory(iDomain).x = interp1(tData,trajectory(iDomain).x,t);
                        trajectory(iDomain).t = t';
                        trajectory(iDomain).z = repmat(trajectory(iDomain).z(1,:)',[1,length(t)])';
                    else % if |tDomain| < dt 
                        iDomain_del = [iDomain_del, iDomain]; % mark row to be deleted
                    end
                end
                trajectory(iDomain_del) = []; % delete rows with timeframe shorter than dt
            end
            
            if ~isequal(obj.systemRef.flipOperator,eye(obj.systemRef.nState))
                trajectory = obj.getCompleteTrajectory(trajectory);
            end
            % loop over trajectories and get constraint forces and
            % actuation
            for iDomain = 1:length(trajectory)
                t  = trajectory(iDomain).t;
                x  = trajectory(iDomain).x;
                z  = trajectory(iDomain).z;
                n  = length(t);
                
                trajectory(iDomain).tau    = zeros(n,obj.systemRef.model.nTau);
                trajectory(iDomain).f      = zeros(n,obj.systemRef.nState);
                trajectory(iDomain).lambda = zeros(n,3,obj.systemRef.model.nZ);
                
                for i = 1:n
                    trajectory(iDomain).tau(i,:) = obj.systemRef.controller.inputTau(t(i),x(i,:)',z(i,:)',obj.systemRef.model,T)';
                    [f, lambda] = obj.systemRef.model.FlowMap(t(i),x(i,:)',z(i,:)',T, obj.systemRef.controller, false, false, true);
                    trajectory(iDomain).f(i,:) = f';
                    lCount = 0;
                    for iZ = 1:obj.systemRef.model.nZ
                        if z(i,iZ) 
                            nLambda = obj.systemRef.model.nLambda(iZ);
                            trajectory(iDomain).lambda(i,1:nLambda,iZ) = lambda(lCount+(1:nLambda))';
                            lCount = lCount + nLambda;
                        end
                    end
                end
            end
        end    
    end
    
    methods (Access = public) % shooting methods
       [xT,simData,trajectory] = singleShootingTD(obj,x0,tDomain,Grad1,Grad2,epsFix) 
       [xT,simData,trajectory] = singleShootingED(obj,x0,T,Grad1,Grad2,epsFix) 
    end
    methods (Access = private) % handles for odeSolver
       dX = DomainDynamics(obj,t,X,z,T,Grad1,Grad2,epsFix,integrandHandle)
       [integrand,tau,tau_t,tau_x,tau_xi,tau_T] = computeIntegrands(obj,t,x,z,T,nInt,integrandHandle)
       trajectory = getCompleteTrajectory(obj,trajectory,tDomain)
    end
end