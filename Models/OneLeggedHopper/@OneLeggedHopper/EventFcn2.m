function [fun, grad_t, grad_x, grad_xi, grad_eps, grad_T] = EventFcn2(obj, t, x, T, controller, Grad1, Grad2, epsFix)
%EventFcn2 lift-off event of 2D hopper
%   Event Lift-Off is defined kinetically by the scaled normal force
%   The gradients of fun w.r.t. x, xi and eps are also provided if
%   Grad1=true

z       = [1;0]; % in stance before lift-off
epsilon = obj.getEpsilon();
k_l     = obj.parameters.param_k_l(epsilon);
l_0     = obj.parameters.l_0;

[~,lambda,~, lambdaGrad] = obj.FlowMap(t, x, z, T, controller, Grad1, Grad2, epsFix);

F_l = k_l*(l_0-x(4)); % joint force

fun = epsilon*lambda(2) + (1-epsilon)*F_l*cos(x(3));

if ~Grad1 
    grad_x   = 0;
    grad_xi  = 0;
    grad_eps = 0;
    grad_t   = 0;
    grad_T   = 0;
else
    F_l_x   = [0,0,0,-k_l, 0,0,0,0];
    
    grad_x  = epsilon*lambdaGrad.lambda_x(2,:)+(1-epsilon)*F_l_x*cos(x(3))...
              -(1-epsilon)*F_l*sin(x(3))*[0 0 1 zeros(1,5)];
    grad_xi = epsilon*lambdaGrad.lambda_xi(2,:);
    grad_eps = [];
    grad_t   = 0;
    grad_T   = 0;
end
end

