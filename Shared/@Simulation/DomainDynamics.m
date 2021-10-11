function dX = DomainDynamics(obj,t,X,z,T,Grad1,Grad2,epsFix,integrandHandle)
    % set up
    nState     = obj.systemRef.nState;
    nCtrl      = obj.systemRef.nCtrl;
    nInt       = length(integrandHandle);  
    if epsFix
        nParam = nCtrl; % don't include nInt here
        nEps   = 0;
    else
        nParam = nCtrl + 1; % including epsilon 
        nEps   = 1;
    end
    if obj.systemRef.controller.timeBased
        nT     = 2;
        nParam = nParam + 2; % including time and period T
    else
        nT     = 0;
    end
    
    x             = X(1:nState);
    [f, ~, fGrad] = obj.systemRef.model.FlowMap(t, x, z, T, obj.systemRef.controller, Grad1, Grad2, epsFix);
    % compute integrands
    [integrand,tau,tau_t,tau_x,tau_xi,tau_T] = obj.computeIntegrands(t,x,z,T,nInt,integrandHandle);
 
    if Grad1  
        shiftIDX = nState+nInt; % shifts index in X
        % linear matrix ODE
        Phi      = reshape(X(shiftIDX+(1:nState^2)),[nState,nState]); % reshape column-wise
        dPhi     = fGrad.f_x*Phi;
        shiftIDX = shiftIDX+nState^2; % update: shifts index in X
             
        % in control parameter
        x_xi     = reshape(X(shiftIDX+(1:nState*nCtrl)),[nState,nCtrl]); % reshape column-wise
        dx_xi    = fGrad.f_x*x_xi + fGrad.f_xi;
        shiftIDX = shiftIDX+nState*nCtrl; % update: shifts index in X
        
        if nT==2
            x_tT     = reshape(X(shiftIDX+(1:nState*nT)),[nState,nT]); % reshape column-wise
            dx_tT    = fGrad.f_x*x_tT + [fGrad.f_t,fGrad.f_T];
            shiftIDX = shiftIDX+nState*nT; % update: shifts index in X
        else
            dx_tT     = [];
        end
        
        if ~epsFix
            x_eps   = reshape(X(shiftIDX+(1:nState*nEps)),[nState,nEps]); % reshape column-wise
            dx_eps  = fGrad.f_x*x_eps + fGrad.f_eps;
        else
            dx_eps  = [];
        end
        dx_p = [dx_xi,dx_tT,dx_eps];
        
        % integrate partial derivatives of user-defined integrands
        dInt_x   = zeros(nInt*nState,1);
        dInt_xi  = zeros(nInt*nCtrl,1);
        dInt_tT  = zeros(nInt*nT,1);
        dInt_eps = zeros(nInt*nEps,1);
        for iInt = 1:nInt
            int_tau       = integrandHandle{iInt}.intF_tau(x,z,tau,t,obj.systemRef);
            int_x         = integrandHandle{iInt}.intF_x(x,z,tau,t,obj.systemRef)+int_tau*tau_x;
            idx_x         = (iInt-1)*nState+(1:nState); % index for deriv. w.r.t. state x
            dInt_x(idx_x) = (int_x*Phi)';
            
            idx_xi          = (iInt-1)*nCtrl + (1:nCtrl); % index for deriv. w.r.t. parameter xi
            int_xi          = int_tau*tau_xi;
            dInt_xi(idx_xi) = (int_x*x_xi + int_xi)';
            
            if nT==2
                int_t         = int_tau*tau_t;
                dInt_tT(iInt) = (int_x*x_tT(:,1) + int_t)';  
                int_T         = int_tau*tau_T;
                dInt_tT(nInt+iInt) = (int_x*x_tT(:,2) + int_T)';  
            end
            
            if ~epsFix
                int_eps        = integrandHandle{iInt}.intF_eps(x,z,tau,t,obj.systemRef);
                dInt_eps(iInt) = (int_x*x_eps + int_eps)';
            end
        end
        dInt = [dInt_x;dInt_xi;dInt_tT;dInt_eps];
        
        
        dX = [f;...
              integrand;...
              reshape(dPhi,[nState^2,1]);... % reshape column-wise
              reshape(dx_p,[nState*nParam,1]);... % reshape column-wise
              dInt];
    else
        dX = [f;integrand];
    end
end