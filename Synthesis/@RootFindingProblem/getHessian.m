function Lxx = getHessian(obj,rootfun,x,lambda)
%getHessian returns second derivative of Lagrangian if
%useLagrangeMultiplier, else compute nullspace derivative

if obj.options.Grad2
    error('Analytic Hessian must be implemented first!')
else    
    Lxx = HessianFD(obj,rootfun,x,lambda);
end

end

function Lxx = HessianFD(obj,rootfun,x,lambda)

stepSize = obj.options.FiniteDifferenceStepSize;

% deactivate optimal control problem for gradient computation
obj.setOptimalControl(false);
obj.updateNumDecVar();

% get f(x)
[~,Lxx] = getJacobianFD(@(x) getNecCon(obj,rootfun,x,lambda),x,stepSize);

%% reset
% reset controller
obj.controller.setFreeParameters(x(obj.idx.Ctrl));
% reset epsilon
if ~isempty(obj.idx.Eps)
    obj.model.setEpsilon(x(obj.idx.Eps));
end
obj.setOptimalControl(true);
obj.updateNumDecVar();

end

function necCon = getNecCon(obj,rootfun,x,lambda)
% get necessary condition

if obj.options.Grad1
       [~,jacobian,objective] = rootfun(x);
       dcost = objective.grad1;
else
       stepSize     = obj.options.FiniteDifferenceStepSize;
       [~,jacobian] = getJacobianFD(rootfun,x,stepSize);
       
       [~,dcost]    = getJacobianFD(@(x) getCostValue(rootfun,x),x,stepSize);     
end

       function cost = getCostValue(rootfun,x)
            [~,~,objec] = rootfun(x);
            cost        = objec.cost;
       end

   if ~isempty(lambda)
        necCon = dcost(1:obj.idx.Ctrl(end))' + jacobian(:,1:obj.idx.Ctrl(end))'*lambda;
   else
        necCon = getTangent(jacobian(:,1:obj.idx.Ctrl(end)))'*dcost(1:obj.idx.Ctrl(end))';
   end
end
