function trajectories = getContinuationTrajectories(LC,type,conParRange,conParSize,tSize,direction)
    if nargin < 6
        direction = 1;
    end
    
    if ischar(type)
        typeChar = type;
    else
        % type is initial state
        typeChar = ['x',num2str(type)];
    end
    
    % make sure that LC is starting Limit Cycle
    numLC  = LC.treeDepth() +1;
    nState = LC.nState;

    nDomain    = length(tSize);
    conParNew  = linspace(conParRange(1),conParRange(2),conParSize);
    
    trajectories = struct('t',cell(conParSize,nDomain),'x',cell(conParSize,nDomain),typeChar,cell(conParSize,nDomain));
    dataOld      = cell(1,nDomain);
    % allocate space
    for iDomain = 1:nDomain
        dataOld{iDomain} = zeros(numLC,(nState+1)*tSize(iDomain));
    end
    conParOld = zeros(numLC,1);
    
    if direction == 1
        rangeLC = (1:numLC);
    else
       rangeLC = (numLC:-1:1);
    end
    for iLC = rangeLC
        traj = LC.getTrajectory();
        trajData = interpT(traj,tSize);
        for iDomain = 1:nDomain
            dataOld{iDomain}(iLC,:) = trajData{iDomain};
        end
        switch type
            case 'epsilon'
                conParOld(iLC)   = LC.Sol.epsilon;
            case num2cell(1:nState)
                conParOld(iLC)   = LC.Sol.x0(type);
            otherwise
                conParOld(iLC)   = LC.Sol.param.(type);
        end       
        %update
        if direction == 1
            LC = LC.Next;
        else
            LC = LC.Prev;
        end
    end
    
    if conParOld(1)~=conParRange(1)
       warning('The starting Limit Cycle might not be correct!') 
    end

    % interpolate
    for iDomain = 1:nDomain
        traj = interp1(conParOld,dataOld{iDomain},conParNew,'spline');
        for iConPar = 1:conParSize
           N = tSize(iDomain);
           trajectories(iConPar,iDomain).t = traj(iConPar,1:N)';
           trajectories(iConPar,iDomain).x = reshape(traj(iConPar,N+1:end),[N,nState]);
           trajectories(iConPar,iDomain).(typeChar) = conParNew(iConPar);
        end
    end

    function trajDataCell = interpT(trajDataIN,tSize)
        nDom = length(tSize);
        trajDataCell = cell(1,nDom);
        for iDom = 1:nDom
            tIN = trajDataIN(iDom).t;
            xIN  = trajDataIN(iDom).x;
            tOUT = linspace(tIN(1),tIN(end),tSize(iDom))';
            xOUT = interp1(tIN,xIN,tOUT,'spline');
            trajDataCell{iDom} = [tOUT',xOUT(:)'];
        end
    end
end