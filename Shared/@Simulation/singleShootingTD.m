function [xT,simData,trajectory] = singleShootingTD(obj,x0,tDomain,Grad1,Grad2,epsFix) 
%singleShootingTD Time-Driven single shooting
%   Detailed explanation goes here
%   !!! the computation of second order gradients (Grad2=true) has not been
%   implemented yet !!!

%% SET UP FOR LOOP OVER DOMAINS
% get dimensions of problem
nTime   = obj.systemRef.nTime;  % number of domain times
nState  = obj.systemRef.nState; % number of states
nCtrl   = obj.systemRef.nCtrl;  % number of control parameters
nInt    = obj.systemRef.nInt;   % number of additional integral-functionals  

% get model, controller, sequence for evaluation
model      = obj.systemRef.model;
controller = obj.systemRef.controller;
sequence   = obj.systemRef.sequence;

if nInt>0
   fncs = obj.systemRef.functionals;
   idx  = cellfun(@(x) x.int,fncs);
   integrandHandle = fncs(idx);
else
   integrandHandle = []; 
end

if epsFix
    nEps   = 0;
    nParam = nCtrl;
else
    nEps   = 1;
    nParam = nCtrl + nEps; % including epsilon 
end

if controller.timeBased
    nT     = 2;
    nParam = nParam + nT; % including time and period T
else
    nT     = 0;
end

% allocate space for flowMapData and jumpMapData
simData = allocateSimData();
% allocate space for trajectory struct
trajectory = struct('t',cell(1,nTime),'x',cell(1,nTime),'z',cell(1,nTime));

xPlus   = x0;
intPlus = zeros(nInt,1);
z       = sequence.FlowMapOrder(:,1);
T       = sum(tDomain);

% loop over domains
for iDomain = 1:nTime
    %% UPDATE XS
    if Grad1 && Grad2
       error('needs to be implemented first!')
    elseif Grad1 && ~Grad2
       XS = [xPlus;...
             intPlus;...
             reshape(eye(nState),[nState^2,1]);... % reshape column-wise
             zeros(nState*nParam,1);...
             zeros((nState+nParam)*nInt,1)];
    else
       XS = [xPlus;intPlus];
    end
    
    %% FLOWMAP
    % get active contacts in FlowMap
    z      = sequence.FlowMapOrder(:,iDomain); 
    odeFun = @(t,X) obj.DomainDynamics(t,X,z,T,Grad1,Grad2,epsFix,integrandHandle);
    % simulate
    t0     = sum(tDomain(1:iDomain-1));
    [t,X]  = obj.odeSolver(odeFun,t0+[0,tDomain(iDomain)],XS,obj.odeOptions);
    xMinus = X(end,1:nState)';
    
    %% JUMPMAP
    % get activated contacts of next domain
    z_prev  = z;
    z_next  = sequence.JumpMapOrder(:,iDomain);
    [xPlus, ~, jumpGrad] = model.JumpMap(xMinus, z_prev, z_next, Grad1, Grad2, epsFix);
    % no jumps in user-defined integrals
    intPlus = X(end,nState+(1:nInt))';
    
    %% WRITE DATA
    % save data
    writeSimData();

    trajectory(iDomain).t = t;
    trajectory(iDomain).x = X(:,1:nState);
    trajectory(iDomain).z = repmat(z,[1,length(t)])';
end
% apply linear flip Operator to exploit symmetry properties
xT      = obj.systemRef.flipOperator*xPlus;


%% nested functions - allocateMapData & writeMapData

    function simData = allocateSimData()
        % integralT   = intT
        % flowMapData = fMD
        % jumpMapData = jMD
        fMD.x_S         = zeros(nState,nTime); % states at start of flowMaps
        fMD.x_E         = zeros(nState,nTime); % states at end of flowMaps
        fMD.integral_S  = zeros(nInt,nTime);   % integrals at start of flowMaps
        fMD.integral_E  = zeros(nInt,nTime);   % integrals at end of flowMaps
        
        if ~Grad1
            % the following is only important to save gradient information
            % and thus can be skipped
            simData.intT        = zeros(nInt,1);
            simData.flowMapData = fMD;
            return
        end
        
        fMD.f_S         = zeros(nState,nTime); % vector-fields at start of flowMaps
        fMD.f_E         = zeros(nState,nTime); % vector-fields at end of flowMaps
        fMD.integrand_S = zeros(nInt,nTime);   % integrands at start of flowMaps
        fMD.integrand_E = zeros(nInt,nTime);   % integrands at end of flowMaps
        % gradients at end of flowMaps w.r.t x_S and p (parameters)
        fMD.xGrad1_E       = zeros(nState,nState,nTime); % = Phi_i(t_i)
        fMD.xiGrad1_E      = zeros(nState,nCtrl,nTime);
        if nT==2
            fMD.tGrad1_E       = zeros(nState,nTime);
            fMD.TGrad1_E       = zeros(nState,nTime);
            fMD.int.tGrad1_E   = zeros(nInt,nTime);
            fMD.int.TGrad1_E   = zeros(nInt,nTime);
        end
        fMD.epsGrad1_E     = zeros(nState,nEps,nTime);
        fMD.int.xGrad1_E   = zeros(nInt,nState,nTime);
        fMD.int.xiGrad1_E  = zeros(nInt,nCtrl,nTime);
        fMD.int.epsGrad1_E = zeros(nInt,nEps,nTime);
        % fMD.xGrad2_E = zeros(nState,nState,nState,nTime);
        % fMD.pGrad2_E = zeros(nState,nParam,nParam,nTime);
        % gradients of reset Maps w.r.t x and p (parameters)
        jMD.xGrad1   = zeros(nState,nState,nTime); % = D_ij
        % control parameters do not affect jump maps, thus:
        jMD.epsGrad1 = zeros(nState,nEps,nTime); 
        % jumpMapData.xGrad2 = zeros(nState,nState,nState,nTime);
        % jumpMapData.pGrad2 = zeros(nState,nParam,nParam,nTime); 
        
        simData.intT        = zeros(nInt,1);
        simData.flowMapData = fMD;
        simData.jumpMapData = jMD;
    end

    function writeSimData()
        XE = X(end,:)';
        tS = t(1);          % start-time of integration
        tE = t(end);        % end-time of integration
        xS = XS(1:nState);  % start state at tS
        xE = XE(1:nState);  % end state at tE
        
        fMD = simData.flowMapData;
        
        if Grad1
            jMD = simData.jumpMapData;
            % avoid clutter data
            clear simData;
        end

        fMD.x_S(:,iDomain)        = xS;
        fMD.x_E(:,iDomain)        = xE;
        fMD.integral_S(:,iDomain) = XS(nState + (1:nInt)); % user-defined integrals
        fMD.integral_E(:,iDomain) = XE(nState + (1:nInt)); % user-defined integrals
        % update terminal integral value
        simData.intT = fMD.integral_E(:,iDomain);
        
        if ~Grad1
            % the following is only important to save gradient information
            % and thus can be skipped
            simData.flowMapData = fMD;
            return
        end
        
        fMD.f_S(:,iDomain)         = model.FlowMap(tS,xS,z,T,controller,false,false,true);
        fMD.f_E(:,iDomain)         = model.FlowMap(tE,xE,z,T,controller,false,false,true);
        fMD.integrand_S(:,iDomain) = obj.computeIntegrands(tS,xS,z,T,nInt,integrandHandle);
        fMD.integrand_E(:,iDomain) = obj.computeIntegrands(tE,xE,z,T,nInt,integrandHandle);
        if Grad1
            shiftIDX = nState+nInt; % shifts index in X
            
            n1 = nState; % number of rows of Jacobian
            n2 = nState; % number of columns of Jacobian
            % = Phi_i(t_i)
            fMD.xGrad1_E(:,:,iDomain)      = reshape(XE(shiftIDX +(1:n1*n2)),[n1,n2]); % reshape column-wise
            shiftIDX = shiftIDX + n1*n2; % update
            
            n1 = nState;
            n2 = nCtrl;
            fMD.xiGrad1_E(:,:,iDomain)     = reshape(XE(shiftIDX +(1:n1*n2)),[n1,n2]); % reshape column-wise 
            shiftIDX = shiftIDX + n1*n2; % update
            
            if nT==2
                n1 = nState;
                n2 = nT;
                tTGrad1_E = reshape(XE(shiftIDX +(1:n1*n2)),[n1,n2]); % reshape column-wise 
                fMD.tGrad1_E(:,iDomain) = tTGrad1_E(:,1);
                fMD.TGrad1_E(:,iDomain) = tTGrad1_E(:,2);
                shiftIDX = shiftIDX + n1*n2; % update
            end
            
            if ~epsFix
                n1 = nState;
                n2 = nEps;
                fMD.epsGrad1_E(:,:,iDomain) = reshape(XE(shiftIDX +(1:n1*n2)),[n1,n2]); % reshape column-wise 
                shiftIDX = shiftIDX + n1*n2; % update
            end
            
            if nInt>0
                n1 = nInt;
                n2 = nState;
                fMD.int.xGrad1_E(:,:,iDomain)  = reshape(XE(shiftIDX +(1:n1*n2)),[n2,n1])'; % reshape row-wise
                shiftIDX = shiftIDX + n1*n2; % update

                n1 = nInt;
                n2 = nCtrl;
                fMD.int.xiGrad1_E(:,:,iDomain) = reshape(XE(shiftIDX +(1:n1*n2)),[n2,n1])'; % reshape row-wise
                shiftIDX = shiftIDX + n1*n2; % update
                
                if nT==2
                    n1 = nInt;
                    n2 = 1;
                    fMD.int.tGrad1_E(:,iDomain) = reshape(XE(shiftIDX +(1:n1*n2)),[n2,n1])'; % reshape row-wise
                    shiftIDX = shiftIDX + n1*n2; % update  

                    fMD.int.TGrad1_E(:,iDomain) = reshape(XE(shiftIDX +(1:n1*n2)),[n2,n1])'; % reshape row-wise
                    shiftIDX = shiftIDX + n1*n2; % update  
                end
                
                if ~epsFix
                    n1 = nInt;
                    n2 = nEps;
                    fMD.int.epsGrad1_E(:,:,iDomain) = reshape(XE(shiftIDX +(1:n1*n2)),[n2,n1])'; % reshape row-wise
                    %shiftIDX = shiftIDX + n1*n2; % update
                end
            end
            
            % jump map data
            jMD.xGrad1(:,:,iDomain)= jumpGrad.xP_x;
            if ~epsFix
                jMD.epsGrad1(:,:,iDomain) = jumpGrad.xP_eps;
            end
        end
        
        simData.flowMapData = fMD;
        simData.jumpMapData = jMD;
    end

end