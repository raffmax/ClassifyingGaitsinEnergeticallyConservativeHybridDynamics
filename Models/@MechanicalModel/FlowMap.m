function [f, lambda, fGrad, lambdaGrad] = FlowMap(obj, t, x, z, T, controller, Grad1, Grad2, epsFix)
%FlowMap General implementation of the vector-field computation of a 
%mechanical model. 
%   FlowMap computes the right-hand side of the first-order ode at time-
%   instance t.

q       = x(1:obj.nQ);
dq      = x(obj.nQ+(1:obj.nQ));   
epsilon = obj.getEpsilon();

% -- function calls --
% - f(x, t) components at current instance -
M = obj.M(epsilon, q);
h = obj.h(epsilon, q, dq, z);
B = obj.B(epsilon, q, dq, z);

activeConstraints = find(z); % only evaluate active contacts
if ~isempty(activeConstraints)
    nActive = length(activeConstraints);
    Wdyn = [];
    W    = [];
    Wdot = [];
    for iActive = 1:nActive
        % get active contact Jacobians
        Wdyn  = [Wdyn,obj.Wdyn{activeConstraints(iActive)}(epsilon, q, z)]; %#ok<AGROW>
        W     = [W,obj.W{activeConstraints(iActive)}(epsilon, q, z)]; %#ok<AGROW>
        Wdot  = [Wdot,obj.Wdot{activeConstraints(iActive)}(epsilon, q, dq, z)]; %#ok<AGROW>
    end
    
    activeHomotopyCon = obj.homotopyCon(:) & z(:); % find active homotopy constraints
    activeHomotopyCon = activeHomotopyCon(activeConstraints);
end

% update controller
controller.update(B,M);
[u,   u_t,   u_x,   u_xi,   u_T]   = controller.inputU( t, x, z, obj, T);
[tau, tau_t, tau_x, tau_xi, tau_T] = controller.inputTau( t, x, z, obj, T);

% -- dynamics --
% - compute contact forces -
% unconstrained flow: M*q''= h+B*tau
ddqf = M\(h + B*tau);

if ~isempty(activeConstraints)
    % Delassus matrix
    QG = M\Wdyn;     % QG corresponds to M^(-1)*Wdyn 
    G  = W'*QG;

    % compute contact force
    lambda = -G\(Wdot'*dq + W'*ddqf);
    
    % homotopy constraint
    lambda(activeHomotopyCon) = (1-epsilon)*lambda(activeHomotopyCon);

    % contact dynamic
    % constrained flow: M*q''= Wdyn*lambda
    ddqc = M\(Wdyn*lambda);
else
    ddqc = zeros(size(ddqf));
    lambda = 0;
end

% combine results to vector field f
ddq   = ddqf + ddqc;
f     = [dq; ddq]+u;

if Grad1
    % --- First Order Sensitivity Analysis ---
    % - derivatives -
    h_x        = obj.Grad1.h_x(epsilon, x, z);     
    B_x_tau    = obj.Grad1.B_x(epsilon, x, z, tau);
    M_x_ddq    = obj.Grad1.M_x(epsilon, x, ddq);


    if ~isempty(activeConstraints)
        Wdyn_x_lambda = zeros(obj.nQ,2*obj.nQ);
        W_x_T_ddq     = [];
        Wdot_x_T_dq   = [];
        lambCount      = 0;
        for iActive = 1:nActive
            nLambda = obj.nLambda(activeConstraints(iActive));
            % get active contact Jacobians
            Wdyn_x_lambda = Wdyn_x_lambda+obj.Grad1.Wdyn_x{activeConstraints(iActive)}(epsilon, x, z, lambda(lambCount+1:nLambda));
            W_x_T_ddq     = [W_x_T_ddq;obj.Grad1.W_x_T{activeConstraints(iActive)}(epsilon, x, z, ddq)]; %#ok<AGROW>
            Wdot_x_T_dq   = [Wdot_x_T_dq;obj.Grad1.Wdot_x_T{activeConstraints(iActive)}(epsilon, x, z, dq)]; %#ok<AGROW>
            lambCount     = lambCount + nLambda;
        end
        Wdot__dq_x  = Wdot'*[zeros(obj.nQ), eye(obj.nQ)];
    else
        Wdyn_x_lambda = 0;
    end

    % compute d/dx Btau
    Btau_x = B_x_tau + B*tau_x;

    dyn_x = -M_x_ddq +  h_x + Btau_x + Wdyn_x_lambda; % = Mq''-Wdyn*lambda_x = dyn_x

    if ~isempty(activeConstraints)
        Qd       = M'\W;
        lambda_x = -G\(Wdot_x_T_dq + Wdot__dq_x+ W_x_T_ddq + Qd'*dyn_x);
        % homotopy constraint
        lambda_x(activeHomotopyCon,:) = (1-epsilon)*lambda_x(activeHomotopyCon,:);
        
        ddq_x    = M\(dyn_x + Wdyn*lambda_x);
    else
        ddq_x    = M\dyn_x;
        lambda_x = [];
    end

    % combine results
    f_x = [[zeros(obj.nQ), eye(obj.nQ)]; ...
           ddq_x] + u_x;

    % compute Grad1 w.r.t. to the control parameters
    ddqf_xi = M\(B*tau_xi);
    if ~isempty(activeConstraints)
        lambda_xi = -G\(W'*ddqf_xi);
        % homotopy constraint
        lambda_xi(activeHomotopyCon,:) = (1-epsilon)*lambda_xi(activeHomotopyCon,:);
        
        ddqc_xi   = M\(Wdyn*lambda_xi);
        ddq_xi    = ddqf_xi + ddqc_xi;
    else
        ddq_xi    = ddqf_xi;
        lambda_xi = [];
    end
    f_xi = [zeros(size(ddq_xi)); ddq_xi]+u_xi;
    
    % compute Grad1 w.r.t. to period T
    if controller.timeBased
        ddqf_t = M\(B*tau_t);
        if ~isempty(activeConstraints)
            lambda_t = -G\(W'*ddqf_t);
            % homotopy constraint
            lambda_t(activeHomotopyCon) = (1-epsilon)*lambda_t(activeHomotopyCon);
            
            ddqc_t   = M\(Wdyn*lambda_t);
            ddq_t    = ddqf_t + ddqc_t;
        else
            lambda_t = [];
            ddq_t    = ddqf_t;
        end
        f_t = [zeros(size(ddq_t)); ddq_t]+u_t;
        
        ddqf_T = M\(B*tau_T);
        if ~isempty(activeConstraints)
            lambda_T = -G\(W'*ddqf_T);
            % homotopy constraint
            lambda_T(activeHomotopyCon) = (1-epsilon)*lambda_T(activeHomotopyCon);
            
            ddqc_T   = M\(Wdyn*lambda_T);
            ddq_T    = ddqf_T + ddqc_T;
        else
            lambda_T = [];
            ddq_T    = ddqf_T;
        end
        f_T = [zeros(size(ddq_T)); ddq_T]+u_T;
    else
        f_t = [];
        f_T = [];
        lambda_t = [];
        lambda_T = [];
    end
            
    %  - compute Grad1 w.r.t. epsilon -
    if ~epsFix
        h_eps     = obj.Grad1.h_eps(epsilon, x, z);     
        Btau_eps  = obj.Grad1.B_eps(epsilon, x, z, tau); % d/d(eps) Btau
        M_eps_ddq = obj.Grad1.M_eps(epsilon, x, ddq);

        if ~isempty(activeConstraints)
            Wdyn_eps_lambda = zeros(obj.nQ,1);
            W_eps_T_ddq     = [];
            Wdot_eps_T_dq   = [];
            lambCount      = 0;
            for iActive = 1:nActive
                nLambda = obj.nLambda(activeConstraints(iActive));
                % get active contact Jacobians
                Wdyn_eps_lambda = Wdyn_eps_lambda+obj.Grad1.Wdyn_eps{activeConstraints(iActive)}(epsilon, x, z, lambda(lambCount+1:nLambda));
                W_eps_T_ddq     = [W_eps_T_ddq;obj.Grad1.W_eps_T{activeConstraints(iActive)}(epsilon, x, z, ddq)]; %#ok<AGROW>
                Wdot_eps_T_dq   = [Wdot_eps_T_dq;obj.Grad1.Wdot_eps_T{activeConstraints(iActive)}(epsilon, x, z, dq)]; %#ok<AGROW>
                lambCount       = lambCount + nLambda;
            end
            % Wdot__dq_eps  = 0
        else
            Wdyn_eps_lambda = 0;
        end

        dyn_eps  = -M_eps_ddq + h_eps + Btau_eps + Wdyn_eps_lambda; % = Mq''-Wdyn*lambda_eps = dyn_eps

        if ~isempty(activeConstraints)
            %Qd = M'\W;
            lambda_eps = -G\(Wdot_eps_T_dq + W_eps_T_ddq + Qd'*dyn_eps);
            
            % homotopy constraint
            add_epsGrad = G\(Wdot'*dq + W'*ddqf);
            lambda_eps(activeHomotopyCon) = (1-epsilon)*lambda_eps(activeHomotopyCon) + add_epsGrad(activeHomotopyCon);
            
            ddq_eps    = M\(dyn_eps + Wdyn*lambda_eps);
        else
            ddq_eps    = M\dyn_eps;
            lambda_eps = [];
        end

        % combine results
        f_eps = [zeros(obj.nQ,1); ddq_eps];
    else
        f_eps      = [];
        lambda_eps = [];
    end
    
    
    fGrad.f_t      = f_t;
    fGrad.f_x      = f_x;
    fGrad.f_xi     = f_xi;
    fGrad.f_T      = f_T;
    fGrad.f_eps    = f_eps;
    fGrad.f_xx     = [];
    fGrad.f_xixi   = [];
    fGrad.f_epseps = [];
    
    lambdaGrad.lambda_t   = lambda_t;
    lambdaGrad.lambda_x   = lambda_x;
    lambdaGrad.lambda_xi  = lambda_xi;
    lambdaGrad.lambda_T   = lambda_T;
    lambdaGrad.lambda_eps = lambda_eps;
else
    fGrad.f_t      = [];
    fGrad.f_x      = [];
    fGrad.f_xi     = [];
    fGrad.f_T      = [];
    fGrad.f_eps    = [];
    fGrad.f_xx     = [];
    fGrad.f_xixi   = [];
    fGrad.f_epseps = [];
    
    lambdaGrad.lambda_t   = [];
    lambdaGrad.lambda_x   = [];
    lambdaGrad.lambda_xi  = [];
    lambdaGrad.lambda_T   = [];
    lambdaGrad.lambda_eps = [];   
end

end