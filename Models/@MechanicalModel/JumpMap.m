function [xPlus, Lambda, jumpGrad] = JumpMap(obj, xMinus, z_prev, z_next, Grad1, Grad2, epsFix)
%JumpMap General implementation of the projection from pre-event-state to 
%post-event-state. 
%   JumpMap computes the discrete mapping at time-instance t_i. 

if z_prev == z_next
    xPlus           = xMinus;
    Lambda          = 0;
    jumpGrad.xP_x   = eye(2*obj.nQ);
    jumpGrad.xP_eps = zeros(2*obj.nQ,1);
    return
end

qM      = xMinus(1:obj.nQ);
dqM     = xMinus(obj.nQ+(1:obj.nQ));
e       = obj.parameters.coefRes;
epsilon = obj.getEpsilon();

% -- function calls --
M = obj.M(epsilon, qM);

activeConstraints = find(z_next); % only evaluate active contacts in next domain
if ~isempty(activeConstraints)
    nActive = length(activeConstraints);
    
    activeHomotopyCon = obj.homotopyCon(:) & z_next(:); % find active homotopy constraints
    activeHomotopyCon = activeHomotopyCon(activeConstraints);
    
    Wdyn = [];
    W    = [];
    for iActive = 1:nActive
        % get active contact Jacobians
        Wdyn  = [Wdyn,obj.Wdyn{activeConstraints(iActive)}(epsilon, qM, z_next)]; %#ok<AGROW>
        W     = [W,obj.W{activeConstraints(iActive)}(epsilon, qM, z_next)]; %#ok<AGROW>
    end

    % Delassus matrix
    QG = M\Wdyn;     % QG corresponds to M^(-1)*Wdyn 
    G  = W'*QG;

    % compute contact force
    Lambda  = -(1+e)*(G\(W'*dqM));
    
    % homotopy constraint
    Lambda(activeHomotopyCon) = (1-epsilon)*Lambda(activeHomotopyCon);
    
    % compute M^(-1)*Wdyn*Lambda
    % i.e. dqPlus-dqMinus
    dqDelta = QG*Lambda;
else % reconstruct lambda
    Lambda  = 0;
    dqDelta = 0*dqM;
end

dqP      = dqM + dqDelta;
xPlus    = [qM; dqP]; % no jump in position!

if Grad1
    if ~isempty(activeConstraints)
        M_x_dqDelta  = obj.Grad1.M_x(epsilon, xMinus, dqDelta);
        %
        Wdyn_x_Lambda = zeros(obj.nQ,2*obj.nQ);
        W_x_T_dqP     = []; % dW/dx' * (dqP + e*dqM)
        lambCount     = 0;
        for iActive = 1:nActive
            nLambda = obj.nLambda(activeConstraints(iActive));
            % get active contact Jacobians
            Wdyn_x_Lambda = Wdyn_x_Lambda+obj.Grad1.Wdyn_x{activeConstraints(iActive)}(epsilon, xMinus, z_next, Lambda(lambCount+1:nLambda));
            W_x_T_dqP     = [W_x_T_dqP;obj.Grad1.W_x_T{activeConstraints(iActive)}(epsilon, xMinus, z_next, dqP + e*dqM)]; %#ok<AGROW>
            lambCount     = lambCount + nLambda;
        end

        Qd     = M'\W;
        jump_x = M_x_dqDelta - Wdyn_x_Lambda;
        dqM_x  = [zeros(obj.nQ), eye(obj.nQ)];

        % compute contact force derivative
        Lambda_x = G\(Qd'*jump_x - (1+e)*W'*dqM_x - W_x_T_dqP);
        
        % homotopy constraint
        Lambda_x(activeHomotopyCon,:) = (1-epsilon)*Lambda_x(activeHomotopyCon,:);

        % combine
        dqP_x    = QG*Lambda_x-M\jump_x + dqM_x;  
        % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if ~epsFix
            dqM_eps        = zeros(obj.nQ,   1);
            M_eps_dqDelta  = obj.Grad1.M_eps(epsilon, xMinus, dqDelta);
            %
            Wdyn_eps_Lambda = zeros(obj.nQ,1);
            W_eps_T_dqP     = []; % dW/depsilon' * (dqP + e*dqM)
            lambCount       = 0;
            for iActive = 1:nActive
                nLambda = obj.nLambda(activeConstraints(iActive));
                % get active contact Jacobians
                Wdyn_eps_Lambda = Wdyn_eps_Lambda+obj.Grad1.Wdyn_eps{activeConstraints(iActive)}(epsilon, xMinus, z_next, Lambda(lambCount+1:nLambda));
                W_eps_T_dqP     = [W_eps_T_dqP;obj.Grad1.W_eps_T{activeConstraints(iActive)}(epsilon, xMinus, z_next, (dqP + e*dqM))]; %#ok<AGROW>
                lambCount       = lambCount + nLambda;
            end

            % Qd = M'\W; % already exists
            jump_eps   = M_eps_dqDelta - Wdyn_eps_Lambda;
            % dqM_eps    = zeros(obj.nQ,2*obj.nQ);
            
            % compute contact force derivative
            Lambda_eps = G\(Qd'*jump_eps - W_eps_T_dqP);
            
            % homotopy constraint
            add_epsGrad = G\(W'*dqM);
            Lambda_eps(activeHomotopyCon) = (1-epsilon)*Lambda_eps(activeHomotopyCon) + add_epsGrad(activeHomotopyCon);

            % combine
            dqP_eps  = QG*Lambda_eps-M\jump_eps + dqM_eps; 
        end
    else
        dqP_x   = [zeros(obj.nQ), eye(obj.nQ)]; % dqM_x  = [zeros(obj.nQ), eye(obj.nQ)];
        dqP_eps = zeros(obj.nQ,1);               % dqM_eps  = zeros(obj.nQ,   1);
    end

    jumpGrad.xP_x = [[eye(obj.nQ), zeros(obj.nQ)]; dqP_x]; % = D_ij
    
    % controller parameter don't affect jump map
    % jumpGrad.xP_xi = zeros(2*obj.nQ, obj.nXi);
    
    if ~epsFix
        jumpGrad.xP_eps = [zeros(obj.nQ,1);dqP_eps];
    else
        jumpGrad.xP_eps = [];
    end
    jumpGrad.xP_xx     = [];
    jumpGrad.xP_epseps = [];
else
    jumpGrad.xP_x      = [];
    jumpGrad.xP_eps    = [];
    jumpGrad.xP_xx     = [];
    jumpGrad.xP_epseps = [];
end
end