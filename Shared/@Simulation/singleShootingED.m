function [xT,simData,trajectory] = singleShootingED(obj,x0,T,Grad1,Grad2,epsFix) 
%singleShootingED Event-Driven single shooting
%   Detailed explanation goes here
%   !!! the computation of second order gradients (Grad2=true) has not been
%   implemented yet !!!

%% SET UP FOR LOOP OVER DOMAINS
% get dimensions of problem
nState  = obj.systemRef.nState; % number of states
nCtrl   = obj.systemRef.nCtrl;  % number of control parameters
nInt    = obj.systemRef.nInt;   % number of additional integral-functionals  

% get model, controller, sequence for evaluation
model      = obj.systemRef.model;
controller = obj.systemRef.controller;
sequence   = obj.systemRef.sequence;

if controller.timeBased
    warning('Period Time needs to be defined as part of the controller!')
end

if nInt>0
   fncs = obj.systemRef.functionals;
   idx  = cellfun(@(x) x.int,fncs);
   integrandHandle = fncs(idx);
else
   integrandHandle = []; 
end

if epsFix
    nEps   = 0;
    nParam = nCtrl;
else
    nEps   = 1;
    nParam = nCtrl + nEps; % including epsilon 
end

% allocate space for trajectory struct
trajectory = struct('t',cell(1,1),'x',cell(1,1),'z',cell(1,1));

simData = [];

xPlus   = x0;
intPlus = zeros(nInt,1);
z       = sequence.FlowMapOrder(:,1);
t       = 0;

%% define events in odeOption
eventIDs = unique(sequence.events);
% get all possible events
options        = obj.odeOptions;

    function [conEvent,isterminal,direction] = events(t,x,z,obj,eventIDs,controller,Grad1,Grad2,epsFix)
        counter  = 0;
        conEvent   = zeros(length(eventIDs),1);
        isterminal = ones(length(eventIDs),1);
        direction  = zeros(length(eventIDs),1);
        for iEvent = eventIDs
            counter = counter+1;
            [conEvent(counter),conData] = obj.systemRef.model.getEventFunctional(t,x,z,[],iEvent,controller, Grad1, Grad2, epsFix);

            if isfield(conData,'direction')
                direction(counter) = conData.direction;
            else
                direction(counter) = 1; % always from - to +
            end
        end
    end

countDomains = 0;
while t(end) < T
    countDomains = countDomains + 1;
    %% UPDATE XS
    if Grad1 && Grad2
       error('needs to be implemented first!')
    elseif Grad1 && ~Grad2
       XS = [xPlus;...
             intPlus;...
             reshape(eye(nState),[nState^2,1]);... % reshape column-wise
             zeros(nState*nParam,1);
             zeros((nState+nParam)*nInt,1)];
    else
       XS = [xPlus;intPlus];
    end
    
    %% FLOWMAP
    odeFun = @(t,X) obj.DomainDynamics(t,X,z,T,Grad1,Grad2,epsFix,integrandHandle);
    if ~isempty(eventIDs)
        options.Events = @(t,x) events(t,x,z,obj,eventIDs,controller,Grad1,Grad2,epsFix);
    end

    % simulate
    if isempty(eventIDs)
        [t,X,] = obj.odeSolver(odeFun,[t(end),T],XS,options);
        ie = 0;
    else
        [t,X,~,~,ie] = obj.odeSolver(odeFun,[t(end),T],XS,options);
    end
    xMinus = X(end,1:nState)';
    
    if ~isempty(ie)
        if isempty(eventIDs)
            eData.z_next = z;
        else
            [~,eData] = model.getEventFunctional(t(end),xMinus,z,[],eventIDs(ie(1)),controller, Grad1, Grad2, epsFix);
        end
        %% JUMPMAP
        % get activated contacts of next domain
        z_prev  = z;
        z_next  = eData.z_next;
        xPlus   = model.JumpMap(xMinus, z_prev, z_next, Grad1, Grad2, epsFix);
        % no jumps in user-defined integrals
        intPlus = X(end,nState+(1:nInt))';

    else
        xPlus  = xMinus;
        z_next = z;
    end
    
    trajectory(countDomains).t = t;
    trajectory(countDomains).x = X(:,1:nState);
    trajectory(countDomains).z = repmat(z,[1,length(t)])';
    
    % update discrete state
    z = z_next;
end
% apply linear flip Operator to exploit symmetry properties
xT      = obj.systemRef.flipOperator*xPlus;
end