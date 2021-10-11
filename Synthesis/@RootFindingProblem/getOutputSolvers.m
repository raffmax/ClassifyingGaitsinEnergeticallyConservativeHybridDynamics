function [sol,output] = getOutputSolvers(obj,init,decVar,output,jacobian,exitflag)
%getOutputSolvers 

% reassamble x0
x0                  = init.x0;
x0(obj.freeStates)  = decVar(obj.idx.FreeState);


tDomain             = decVar(obj.idx.Time);
T                   = sum(tDomain); % Period Time
xi                  = obj.controller.getFreeParameters();
epsilon             = obj.model.getEpsilon();

sol                        = Solution(x0,tDomain,xi,epsilon,init);
sol.rfpData.decVar         = decVar;
sol.rfpData.periodicStates = obj.periodicStates;

if ~isequal(obj.flipOperator, eye(obj.nState))
    sol.rfpData.flipOperator = obj.flipOperator;
else
    sol.rfpData.flipOperator = [];
end

if obj.fixedParameter
    field             = obj.functionals{1}.c;
    sol.param.(field) = init.param.(field);
end

if isfield(output,'t')
    sol.rfpData.jacobianAug = [jacobian;output.t'];
end

obj.fixEpsilon(true);
[~,~,objective] = rootFunctionTDsingle(obj,init,decVar,false);

if ~isempty(objective)
    if isfield(objective,'cost')
        cost = objective.cost;
    else
        cost = [];
    end
    if isfield(objective,'cost')
        costGrad1 = objective.grad1(1:obj.idx.Ctrl(end))';
    else
        costGrad1 = [];
    end
else
    cost      = [];
    costGrad1 = [];
end
sol.rfpData.cost = cost;
output.cost      = cost;
sol.rfpData.OptimalCtrl.gradientCost = costGrad1;

if obj.optimalControl
    % get size of x = [x0,tDomain,xi]
    nX   = obj.idx.Ctrl(end); 
    if obj.useLagrangeMultipliers
        nCon        = length(obj.idx.Multiplier);
        multiplier  = decVar(obj.idx.Multiplier); 
        jacobianCon = jacobian(nX+(1:nCon),1:nX);
        output.tCon = null(jacobianCon); %null space
        output.secondDerivative = output.tCon'*jacobian(1:nX,1:nX)*output.tCon;
    else
        nCon        = obj.nCon;
        jacobianCon = jacobian(end-nCon+1:end,1:nX);
        multiplier  = -jacobianCon'\costGrad1;
        output.tCon = getTangent(jacobianCon); %null space
        output.secondDerivative = jacobian(1:end-nCon,1:nX)*output.tCon;
    end
    output.jacobianCon      = jacobianCon;
    output.multiplier       = multiplier;
  
    sol.rfpData.OptimalCtrl.secondDerivative = output.secondDerivative;
    sol.rfpData.OptimalCtrl.jacobianCon      = output.jacobianCon;
    sol.rfpData.OptimalCtrl.tCon             = output.tCon;
    sol.rfpData.multiplier                   = multiplier;
    
else
    output.jacobianCon = jacobian;
    sol.rfpData.multiplier  = decVar(obj.idx.Multiplier); 
end

if exitflag~=1 || ~obj.options.Grad1
    return
end
%% construct Monodromy matrix and get cost
if strcmp(obj.integrationType,'TD') && strcmp(obj.shootingType,'single')
    fun  = @(x0,tDomain) obj.simulation.singleShootingTD(x0,tDomain,true,false,true);
else
    error('Has not been implemented yet!')
end

[~,simData] = fun(x0,tDomain);
% get eventGradients
nEvents         = length(obj.sequence.events);
eventGradients  = zeros(nEvents,obj.nState);
for iEvent = 1:nEvents
    % get x at event
    xe_i = simData.flowMapData.x_E(:,iEvent);
    te_i = sum(tDomain(1:iEvent));
    [~,eventData] = obj.model.getEventFunctional(te_i,xe_i,[],T,obj.sequence.events(iEvent),obj.controller,true,false,true);
    eventGradients(iEvent,:) = eventData.grad_x;
end

sol.dynamics.MonodromyMatrix = obj.simulation.computeMonodromyMatrix(simData,eventGradients);
sol.dynamics.f0              = simData.flowMapData.f_S(:,1);
sol.dynamics.E_x0            = obj.model.Grad1.E_x(epsilon,x0,obj.sequence.FlowMapOrder(:,1));

%% get eventVelocities
sol.rfpData.eventVelocities = diag(eventGradients*simData.flowMapData.f_E);

end
