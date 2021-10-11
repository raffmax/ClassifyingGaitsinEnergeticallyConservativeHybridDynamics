function [sol,fval,exitflag,output,jacobian] = solveNewtonRFP(obj,init,varargin)
%solveNewtonRFP is a wrapper function for NewtonsMethod (Newton's method).
%   solveNewtonRFP constructs the inital decision variable and function 
%   handle to the rootFunction and calls newtonRaphson with this 
%   information.
%   The result is put in the solution struct sol.

if ~isequal(class(init),'Solution')
    error('init must be an object of class Solution!')
end


% set new control parameter in Controller
obj.controller.setFreeParameters(init.xi);
% set new epsilon parameter in model
obj.model.setEpsilon(init.epsilon);

decVarInit = obj.getDecisionVariable(init);
init.rfpData.decVar = decVarInit; % add decision variable to init

if strcmp(obj.integrationType,'TD') && strcmp(obj.shootingType,'single')
    fun = @(x)rootFunctionTDsingle(obj,init,x,obj.optimalControl);
else
    error('Has not been implemented yet!')
end

[decVar,fval,exitflag,output,jacobian] = NewtonsMethod(fun,decVarInit,obj.options);

[sol,output] = obj.getOutputSolvers(init,decVar,output,jacobian,exitflag);

end
