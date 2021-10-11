function decVar = getDecisionVariable(obj,Sol,varargin)
%getDecisionVariable returns the decision variable.
if strcmp(obj.integrationType,'TD') && strcmp(obj.shootingType,'single')
    if obj.fixedEpsilon
        decVar = [Sol.x0(obj.freeStates);Sol.tDomain;Sol.xi];
    else
        decVar = [Sol.x0(obj.freeStates);Sol.tDomain;Sol.xi;Sol.epsilon];
    end
    if obj.optimalControl
        if obj.useLagrangeMultipliers
            if nargin > 2
                % are Lagrange multipliers provided in varargin?
                lambda = varargin{:};
                if length(lambda)~=obj.nCon
                   error('Size of Lagrange multipliers must be equal to number of constraints!') 
                end
            elseif ~isempty(Sol.rfpData.multiplier)
                lambda = Sol.rfpData.multiplier;
            else
                fun = @(x)rootFunctionTDsingle(obj,Sol,x,false);
                % assume optimal point and reconstruct lambda

                [F,jacobian,objective] = fun(decVar);
                % lambda is unique if decVar is optimal point and jacobian has 
                % full rank 
                lambda = -jacobian(:,1:obj.idx.Ctrl(end))'\objective.grad1(1:obj.idx.Ctrl(end))'; 
                if max(abs(F))>obj.options.FunctionTolerance || rank(jacobian)<obj.nCon
                   warning('Reconstructing Lagrange multipliers may not be accurate!') 
                end
            end   
        else
            lambda = [];
        end
        if obj.fixedEpsilon
            decVar = [decVar;lambda];
        else
            decVar = [decVar(1:end-1);lambda;decVar(end)];
        end
    end
elseif strcmp(obj.shootingType,'multiple')
    error('multiple shooting is not implemented yet')
end
end

