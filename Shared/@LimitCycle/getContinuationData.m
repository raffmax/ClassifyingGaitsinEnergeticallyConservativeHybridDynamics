function conData = getContinuationData(LC)
    %getContinuationData
    % collects data from start LC  to the last LC in its tree
    numLC = LC.treeDepth() +1;
    % define struct (allocate space)
    conData.x0      = zeros(numLC,LC.nState);
    conData.tDomain = zeros(numLC,LC.nTime);
    conData.epsilon = zeros(numLC,1);
    conData.xi      = zeros(numLC,LC.nCtrl);
    conData.direction = zeros(numLC, 1);
    if ~isempty(LC.Sol.dynamics)
        conData.MonodromyMatrix = zeros(LC.nState,LC.nState,numLC);
        conData.f0              = zeros(numLC,LC.nState);
    end
    for fn = fieldnames(LC.Sol.rfpData)'
        if ~strcmp('decVar',fn{1}) && ~strcmp('OptimalCtrl',fn{1}) && ~strcmp('jacobianAug',fn{1})
            conData.(fn{1}) = zeros(numLC,length(LC.Sol.rfpData.(fn{1})));
        end 
    end
    if ~isempty(LC.Sol.param)
        for fn = fieldnames(LC.Sol.param)'
           conData.(fn{1}) = zeros(numLC,1);
        end
    end
    conData.eventVelocities = zeros(numLC,length(LC.sequence.events));
    conData.detAugJacobi    = zeros(numLC,1);
    conData.detSimpleJacobi = zeros(numLC,1);
    if ~isempty(conData.multiplier)
        optimProblem = true;
        conData.detConJacobi = zeros(numLC,1);
        conData.detHessian   = zeros(numLC,1);
        conData.suffCondEig  = zeros(numLC,size(LC.Sol.rfpData.OptimalCtrl.tCon,2));
    else
        optimProblem = false;
    end
    for iLC = 1:numLC
        % loop over limit cycles
        if iLC > 1
            LC = LC.getNextLC();
        end
        conData.x0(iLC,:)         = LC.Sol.x0(:)';
        conData.tDomain(iLC,:)    = LC.Sol.tDomain(:)';
        conData.epsilon(iLC)      = LC.Sol.epsilon;
        conData.xi(iLC,:)         = LC.Sol.xi(:)';
        if isfield(LC.Sol.rfpData,'eventVelocities')
            conData.eventVelocities(iLC,:) = LC.Sol.rfpData.eventVelocities;
        end
        conData.detAugJacobi(iLC)      = det(LC.Sol.rfpData.jacobianAug);
        conData.detSimpleJacobi(iLC)   = det(LC.Sol.rfpData.jacobianAug(1:end-1,1:end-1));
        
        if isfield(LC.Sol.rfpData,'direction')
            if ~isempty(LC.Sol.rfpData.direction)
                conData.direction(iLC) = LC.Sol.rfpData.direction;
            end
        end

        
        if isfield(LC.Sol.rfpData,'cost')
            if ~isempty(LC.Sol.rfpData.cost)
                conData.cost(iLC) = LC.Sol.rfpData.cost;
            end
        end
        
        if ~isempty(LC.Sol.param)
            for fn = fieldnames(LC.Sol.param)'
               conData.(fn{1})(iLC) = LC.Sol.param.(fn{1});
            end
        end
        
        if isfield(conData,'MonodromyMatrix')
            conData.MonodromyMatrix(:,:,iLC) = LC.Sol.dynamics.MonodromyMatrix;
            conData.f0(iLC,:) = LC.Sol.dynamics.f0';
        end
        if optimProblem
            conData.multiplier(iLC,:) = LC.Sol.rfpData.multiplier(:)';
            jCon   = LC.Sol.rfpData.OptimalCtrl.jacobianCon;
            sufCon = LC.Sol.rfpData.OptimalCtrl.secondDerivative;
            conData.detConJacobi(iLC)  = det(jCon*jCon');
            conData.detHessian(iLC)    = det(sufCon);
            conData.suffCondEig(iLC,:) = eig(sufCon);
        end
    end
    
    if isfield(conData,'MonodromyMatrix')
        [conData.FloquetVectors,conData.FloquetMultipliers] = eigenshuffle(conData.MonodromyMatrix,conData.f0);
        conData.FloquetMultipliers = conData.FloquetMultipliers';
    end
end