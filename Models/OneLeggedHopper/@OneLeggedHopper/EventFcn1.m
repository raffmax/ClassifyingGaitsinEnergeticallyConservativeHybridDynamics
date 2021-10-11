function [fun, grad_t, grad_x, grad_xi, grad_eps, grad_T] = EventFcn1(~,~,x,~,controller, ~, ~, ~)
%EventFcn1 touch-down event of hopper
%   Event Touch-Down is defined kinematically by fun = y-l*cos(alpha+phi)
%   The gradients of fun w.r.t. x, xi and eps are also provided 
fun      = x(2)-x(5)*cos(x(3)+x(4));
grad_t   = 0;
grad_T   = 0;
grad_x   = [0, 1, x(5)*sin(x(3)+x(4)), x(5)*sin(x(3)+x(4)), -cos(x(3)+x(4)), zeros(1,5)];
grad_xi  = zeros(1,controller.nXi);
grad_eps = 0;
end