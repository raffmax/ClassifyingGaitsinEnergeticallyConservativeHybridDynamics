function nCon = getNumConstraints(obj)
%getNumConstraints counts all constraints
%   Returns the number of constraints defined in RootFindingProblem

% count periodic constraints
nConP = length(obj.periodicStates);
% count events
nConE = length(obj.sequence.events);
% count implicit constraints
nConF = 0;
nFun  = length(obj.functionals);
for iFun = 1:nFun
    if ~obj.functionals{iFun}.cost
        % if not a cost functional, add constraint!
        nConF = nConF + 1;
    end
end

nCon = nConP+nConE+nConF;

% set index/description of equality constraints
obj.idxCon.Periodicity = (1:nConP);
obj.idxCon.Events      = nConP+(1:nConE);
obj.idxCon.Additional  = nConP+nConE+(1:nConF);

if obj.optimalControl %shift index
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
    if obj.useLagrangeMultipliers
        obj.idxCon.Periodicity = obj.idxCon.Periodicity + nDecVar;
        obj.idxCon.Events      = obj.idxCon.Events + nDecVar;
        obj.idxCon.Additional  = obj.idxCon.Additional + nDecVar;
    else
        obj.idxCon.Periodicity = obj.idxCon.Periodicity + nDecVar - nCon;
        obj.idxCon.Events      = obj.idxCon.Events + nDecVar-nCon;
        obj.idxCon.Additional  = obj.idxCon.Additional + nDecVar-nCon;  
    end
end
end