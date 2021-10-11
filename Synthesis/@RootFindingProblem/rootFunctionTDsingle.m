function [fun,jacobian,objective] = rootFunctionTDsingle(obj,init,decVar,computeLxx)
%rootFunctionTDsingle(obj,init,decVar)
%   This function uses time-driven single shooting in simulation to obtain
%   the rootFunction value

if nargin<4
    % don't compute second derivative in optimization problem
    computeLxx = false;
end

% get initial state, inital domain durations and gradient flags 
x0                  = init.x0;
x0(obj.freeStates)  = decVar(obj.idx.FreeState);
z0                  = obj.sequence.FlowMapOrder(:,1);
zT                  = obj.sequence.JumpMapOrder(:,end);
tDomain             = decVar(obj.idx.Time);
T                   = sum(tDomain);
xi                  = decVar(obj.idx.Ctrl);
Grad1               = obj.options.Grad1;
Grad2               = obj.options.Grad2;
epsFix              = obj.fixedEpsilon; % is epsilon fixed or a variable?

% set new control parameters in Controller
obj.controller.setFreeParameters(xi);
% set/get new homotopy parameter in MechanicalModel
if ~obj.fixedEpsilon
    epsilon   = decVar(obj.idx.Eps);
    obj.model.setEpsilon(epsilon);
else
    epsilon   = obj.model.getEpsilon();
end

%% run simulation
[xT,simData,~]    = obj.simulation.singleShootingTD(x0,tDomain,Grad1,Grad2,epsFix);

% define structure of jacobian
if Grad1
    % allocate space for Jacobian
    jacobian = zeros(obj.nCon,obj.nDecVar);
    % compose gradient information
    [xE,intE] = obj.getGradientsFromSimData(simData);
    % compute discrete map to xT
    xT_xE = obj.flipOperator*simData.jumpMapData.xGrad1(:,:,end);
else
    jacobian = [];
end

%% define periodicity constraints
% compute selection matrix in periodicty of x0
selectMatrix = eye(obj.nState);
for index = obj.periodicStatesHomotopy
    selectMatrix(index,index) = epsilon;
end

conPeriod = xT(obj.periodicStates)-selectMatrix(obj.periodicStates,:)*x0;

% add information to jacobian
add2JacobianPeriodicity();

%% define event constraints
nEvents   = length(obj.idxCon.Events);
conEvent  = zeros(nEvents,1);
eventData = cell(nEvents,1);
for iEvent = 1:nEvents
    % get x at event
    xe_i = simData.flowMapData.x_E(:,iEvent);
    te_i = sum(tDomain(1:iEvent));
    [conEvent(iEvent),eventData{iEvent}] = obj.model.getEventFunctional(te_i,xe_i,[],T,obj.sequence.events(iEvent),obj.controller,Grad1,Grad2,epsFix);
    add2JacobianEvent(iEvent);
end
%% define functional constraints
nFun       = length(obj.idxCon.Additional);
conFun     = zeros(nFun,1);
fncs       = obj.functionals;

if ~isempty(fncs)
    idxInt = find(cellfun(@(x) x.int,fncs)); % index of integral-functionals
    
    % loop over all functionals!
    iConFun = 0; % counter of equality constraint functional
    for iFNC = 1:length(fncs)
        if obj.fixedParameter && iFNC==1
            % user-provided fixedParameter c must be in first functional!
            % get constant parameter from Solution-struct init
            c   = init.param.(fncs{iFNC}.c);
        else
            c   = 0;
        end
        
        if ~fncs{iFNC}.cost
            iConFun = iConFun +1; % increase equality constraint counter
            intPos = find(idxInt==iFNC); % find position of user-defined integral
            [conFun(iFNC),tau,tau_t,tau_x,tau_xi,tau_T,intT] = getFunctionalValue(iFNC);
            add2JacobianFunctional(iFNC,iConFun)
        end
    end
end

%% assamble constraints
fun = [conPeriod;...
       conEvent;...
       conFun];
   
%% get objective/cost data
if ~isempty(obj.functionals)
    idxCost = find(cellfun(@(x) x.cost,obj.functionals));
    if isempty(idxCost)
        % there is no objective/cost information for this problem
        objective = [];
    elseif length(idxCost)>1
        error('Multi-Cost-Functionals are not implemented yet!')
    else
        intPos = find(idxInt==idxCost); %find position of user-defined integral
        [objective.cost,tau,tau_t,tau_x,tau_xi,tau_T,intT] = getFunctionalValue(idxCost); 
        if Grad1
            objective.grad1 = getFunctionalGrad1(idxCost);
        end
    end
else
    objective = [];
end

% compute jacobian including Lagrange multipliers
if computeLxx
    h = fun; % equality constraints
    if obj.useLagrangeMultipliers
        % L   = c(X) + lambda'*h(X)
        % fun = [L_X,L_lambda]' = [L_X';h]
        if obj.fixedEpsilon
            idxX     = 1:obj.nDecVar-obj.nCon;
            X        = decVar(idxX);
            Lambda   = decVar(obj.idx.Multiplier);
            h_X      = jacobian(:,idxX);
            L_XX     = obj.getHessian(@(X)rootFunctionTDsingle(obj,init,X),X,Lambda);
            L_X      = objective.grad1(idxX)'+h_X'*Lambda;
            fun      = [L_X;h];
            jacobian = [L_XX, h_X'; h_X, zeros(obj.nCon)];
        else
            idxX     = [(1:obj.nDecVar-obj.nCon-1),obj.nDecVar];
            X        = decVar(idxX);
            Lambda   = decVar(obj.idx.Multiplier);
            h_X      = jacobian(:,idxX);
            L_XX     = obj.getHessian(@(X)rootFunctionTDsingle(obj,init,X),X,Lambda);
            L        = objective.grad1(idxX)'+h_X'*Lambda;
            % separate epsilon from X
            fun      = [L(1:end-1);h];
            L_xx     = L_XX(:,1:end-1);
            L_xeps   = L_XX(:,end);
            h_x      = jacobian(:,1:obj.nDecVar-obj.nCon-1);
            h_eps    = jacobian(:,end);
            jacobian = [L_xx, h_x'          , L_xeps;...
                        h_x ,zeros(obj.nCon), h_eps];         
        end
    else
        % fun = [null(h'(x))*c'(x),h(x)]'
        if obj.fixedEpsilon 
            X        = decVar;
            h_X      = jacobian;
            firstOpt = getTangent(h_X)'*objective.grad1';
            fun      = [firstOpt;h];
            firstOpt_X  = obj.getHessian(@(X)rootFunctionTDsingle(obj,init,X),X,[]); 
            jacobian = [firstOpt_X; h_X];
        else
            X        = decVar;
            h_x      = jacobian(:,1:end-1);
            firstOpt = getTangent(h_x)'*objective.grad1(1:end-1)';
            fun      = [firstOpt;h];
            firstOpt_x = obj.getHessian(@(X)rootFunctionTDsingle(obj,init,X),X,[]); 
            jacobian = [firstOpt_x; jacobian];
        end       
    end
end

   

%% nested functions 
%                   getFunctionalValue, getFunctionalGrad1
%                   add2JacobianPeriodicity, add2JacobianEvent, add2JacobianFunctional    
%
    function [F,tau,tau_t,tau_x,tau_xi,tau_T,intT] = getFunctionalValue(iFNC)
        if ~isempty(intPos)
            intT       = simData.flowMapData.integral_E(intPos,end);
        else
            intT = [];
        end
        
        switch fncs{iFNC}.type
            case 'start'
                [tau, tau_t, tau_x, tau_xi, tau_T] = obj.controller.inputTau(0, x0, z0, obj.model, T);
                F                                  = fncs{iFNC}.F(x0,tau,c,obj);
            case 'end'
                [tau, tau_t, tau_x, tau_xi, tau_T] = obj.controller.inputTau(T, xT, zT, obj.model, T);
                F                                  = fncs{iFNC}.F(xT,intT,tau,T,c,obj);  
        end
    end

    function grad = getFunctionalGrad1(iFNC)
        %allocate gradient
        grad = zeros(1,size(jacobian,2));
        functional = fncs{iFNC};
        switch functional.type
            case 'start'
                F_tau = functional.F_tau(x0,tau,c,obj);
                F_x   = functional.F_x(x0,tau,c,obj) + F_tau*tau_x;
                grad(obj.idx.FreeState) = F_x(obj.freeStates);
                grad(obj.idx.Ctrl)      = F_tau*tau_xi;
                if ~obj.fixedEpsilon
                    % compute F_eps
                    F_eps = functional.F_eps(x0,tau,c,obj);
                    grad(obj.idx.Eps) = F_eps;
                end

            case 'end'
                F_tau = functional.F_tau(xT,intT,tau,T,c,obj);
                F_x   = functional.F_x(xT,intT,tau,T,c,obj) + F_tau*tau_x;
                if isempty(functional.F_int)
                    F_int = [];
                else
                    F_int = functional.F_int(xT,intT,tau,T,c,obj);
                end
                
                F_x0  = F_x*xT_xE*xE.x0(:,:,end);
                if ~isempty(F_int)
                    F_x0 = F_x0 + F_int*intE.x0(intPos,:,end);
                end
                grad(obj.idx.FreeState) = F_x0(obj.freeStates);
                % get time derivatives
                % compute dF/dt_i
                F_T = functional.F_T(xT,intT,tau,T,c,obj);
                if obj.controller.timeBased
                    F_T = F_T + F_tau*(tau_T+tau_t);
                end
                for iTime = 1:obj.nTime
                   F_ti = F_T+F_x*xT_xE*xE.t(:,iTime,end);
                   if ~isempty(F_int)
                       F_ti = F_ti +F_int*intE.t(intPos,iTime,end);
                   end
                   grad(obj.idx.Time(iTime)) = F_ti;
                end
                % compute dF/dxi
                F_xi  = F_tau*tau_xi + F_x*xT_xE*xE.xi(:,:,end);
                if ~isempty(F_int)
                    F_xi = F_xi + F_int*intE.xi(intPos,:,end);
                end
                grad(obj.idx.Ctrl) = F_xi;
                if ~obj.fixedEpsilon
                   % compute dF/deps
                   F_eps = functional.F_eps(xT,intT,tau,T,c,obj)...
                           + F_x*xT_xE*xE.eps(:,:,end)...
                           + F_x*simData.jumpMapData.epsGrad1(:,:,end);
                   if ~isempty(F_int)
                       F_eps = F_eps + F_int*intE.eps(intPos,:,end);
                   end
                   grad(obj.idx.Eps) = F_eps;
                end 
        end
    end

    function add2JacobianPeriodicity()
        if ~Grad1
            return
        end
        
        % get size of periodic states
        nPer = length(obj.periodicStates);
        % compute dxT/dx0
        xT_x0 = xT_xE*xE.x0(:,:,end);
        jacobian(1:nPer,obj.idx.FreeState) = xT_x0(obj.periodicStates,obj.freeStates)-selectMatrix(obj.periodicStates,obj.freeStates);
        % compute dxT/dt_i
        for iTime = 1:obj.nTime
           xT_ti = xT_xE*xE.t(:,iTime,end);
           jacobian(1:nPer,obj.idx.Time(iTime)) = xT_ti(obj.periodicStates);
        end
        % compute dxT/dxi
        xT_xi = xT_xE*xE.xi(:,:,end);
        jacobian(1:nPer,obj.idx.Ctrl) = xT_xi(obj.periodicStates,:);
        if ~obj.fixedEpsilon
           % compute dxT/deps
           xT_eps = xT_xE*xE.eps(:,:,end)...
                    + obj.flipOperator*simData.jumpMapData.epsGrad1(:,:,end);
           
           % d/deps selectMatrix*x0 = (d/deps selectMatrix) * x0
           % = diag(obj.periodicStatesHomotopy) * x0
           % or equivalently -> x0(~obj.periodicStatesHomotopy) = 0
           selMat_epsx0 = x0;
           selMat_epsx0(~ismember((1:obj.nState),obj.periodicStatesHomotopy)) = 0;
           % d xT / d eps - d (selectMatrix*x0) / d eps
           jacobian(1:nPer,obj.idx.Eps) = xT_eps(obj.periodicStates)-selMat_epsx0(obj.periodicStates);
        end
    end

    function add2JacobianEvent(iEvent)
        if ~obj.options.Grad1
            return
        end

        % get size of periodic states
        nPer = length(obj.periodicStates);
        % compute de_i/dx0
        e_x0 = eventData{iEvent}.grad_x*xE.x0(:,:,iEvent);
        idx  = nPer+(1:(obj.nCon-nPer));
        jacobian(idx(iEvent),obj.idx.FreeState) = e_x0(obj.freeStates); 
        % compute de_i/dt_j
        for iTime = 1:obj.nTime
           e_ti = eventData{iEvent}.grad_x*xE.t(:,iTime,iEvent);
           if obj.controller.timeBased
               e_ti = e_ti + eventData{iEvent}.grad_T;
               if iEvent >= iTime
                   e_ti = e_ti + eventData{iEvent}.grad_t;
               end
           end
           jacobian(idx(iEvent),obj.idx.Time(iTime)) = e_ti;
        end
        % compute de_i/dxi
        jacobian(idx(iEvent),obj.idx.Ctrl) = eventData{iEvent}.grad_x*xE.xi(:,:,iEvent) +...
                                             eventData{iEvent}.grad_xi;
        if ~obj.fixedEpsilon
           % compute de_i/deps
           jacobian(idx(iEvent),obj.idx.Eps) = eventData{iEvent}.grad_x*xE.eps(:,:,iEvent) +...
                                               eventData{iEvent}.grad_eps;
        end
    end

    function add2JacobianFunctional(iFNC,iConFun)
        if ~obj.options.Grad1
            return
        end
        
        grad = getFunctionalGrad1(iFNC);
        jacobian(end-nFun+iConFun,:) = grad; % fill jacobian
    end

end