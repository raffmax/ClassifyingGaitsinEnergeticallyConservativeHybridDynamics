function nDecVar = getNumDecVar(obj)
%getNumDecVar counts decision variables
%   returns number of decision variables of RootFindingProblem

nDecVar = 0;
% get number of free states
if strcmp(obj.shootingType,'single')
    nDecVar = nDecVar + length(obj.freeStates);
elseif strcmp(obj.shootingType,'multiple')
    error('multiple shooting is not implemented yet')
end
% get number of domain times
nDecVar = nDecVar + obj.nTime;
% get number of control parameters
nDecVar = nDecVar + obj.nCtrl;
% get number of lagrange multipliers
if obj.optimalControl && obj.useLagrangeMultipliers
	nDecVar = nDecVar + obj.nCon;
end
% get number of variable homotopy parameters
if ~obj.fixedEpsilon  
	nDecVar = nDecVar +1;
end

% set decVar index
obj.idx.FreeState = (1:obj.nFreeState);
obj.idx.Time      = obj.nFreeState+(1:obj.nTime);
obj.idx.Ctrl      = obj.nFreeState+obj.nTime+(1:obj.nCtrl);
if obj.optimalControl && obj.useLagrangeMultipliers
    obj.idx.Multiplier = obj.nFreeState+obj.nTime+obj.nCtrl+(1:obj.nCon);
else
    obj.idx.Multiplier = [];
end

if obj.fixedEpsilon
    obj.idx.Eps   = [];
else
    if obj.optimalControl && obj.useLagrangeMultipliers
        obj.idx.Eps = obj.nFreeState+obj.nTime+obj.nCtrl+obj.nCon + 1;
    else
        obj.idx.Eps = obj.nFreeState+obj.nTime+obj.nCtrl + 1;
    end
end
end

