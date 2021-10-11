function trajectory = getCompleteTrajectory(obj,trajectory)
    nTime = length(trajectory); %obj.systemRef.nTime;
    flipOperator = obj.systemRef.flipOperator;
    % get periodicStates
    if isa(obj.systemRef,'LimitCycle')
        periodicStates = obj.systemRef.Sol.rfpData.periodicStates;
    else
        periodicStates = obj.systemRef.periodicStates;
    end
    stateIdx  = (1:obj.systemRef.nState);
    offsetIdx = stateIdx(~ismember(stateIdx,periodicStates));
    for iDomain = 1:nTime
        trajectory(iDomain+nTime).t = trajectory(iDomain).t-trajectory(iDomain).t(1)+trajectory(nTime+iDomain-1).t(end);
        trajectory(iDomain+nTime).x = trajectory(iDomain).x*flipOperator';
        % add offset
        trajectory(iDomain+nTime).x(:,offsetIdx) = trajectory(iDomain).x(:,offsetIdx)-trajectory(iDomain).x(1,offsetIdx)+ trajectory(iDomain+nTime-1).x(end,offsetIdx);
        if size(trajectory(iDomain).z,2) == 4
            trajectory(iDomain+nTime).z = trajectory(iDomain).z(:,[2,1,4,3]);
        elseif size(trajectory(iDomain).z,2) == 5
            trajectory(iDomain+nTime).z = trajectory(iDomain).z(:,[2,1,4,3,5]);
        else
            trajectory(iDomain+nTime).z = flip(trajectory(iDomain).z,2);
        end
    end
end