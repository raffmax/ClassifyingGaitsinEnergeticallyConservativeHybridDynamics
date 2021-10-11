function monodromy = computeMonodromyMatrix(obj,simData,eventGradients)
%computeMonodromyMatrix 
% i.e. d x(T) / d x(0)
monodromy = eye(obj.systemRef.nState);
nDomain   = obj.systemRef.nTime;

% compute fundamental and saltation matrices
for iDomain = 1:nDomain
    % Fundamental Solution Matrix right before event
    Phi = simData.flowMapData.xGrad1_E(:,:,iDomain)*monodromy;
    
    f1  = simData.flowMapData.f_E(:,iDomain);
    if iDomain <= size(eventGradients,1)
        n   = eventGradients(iDomain,:)';
    else
        n   = f1;
    end
    D   = simData.jumpMapData.xGrad1(:,:,iDomain);
    if iDomain == nDomain
        f2 = simData.flowMapData.f_S(:,1);
    else
        f2 = simData.flowMapData.f_S(:,iDomain+1);
    end
    
    % Saltation Matrix at event
    dt = n'*f1;
    if dt == 0
        S   = D;
    else
        S   = D + (f2-D*f1)*n'/(n'*f1);
    end
    % update
    monodromy = obj.systemRef.flipOperator*S*Phi;
end

end
