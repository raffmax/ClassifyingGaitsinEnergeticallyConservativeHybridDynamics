function [integrand,tau,tau_t,tau_x,tau_xi,tau_T] = computeIntegrands(obj,t,x,z,T,nInt,integrandHandle)
    systemRef = obj.systemRef;
    integrand = zeros(nInt,1);
    [tau, tau_t, tau_x, tau_xi, tau_T] = systemRef.controller.inputTau(t, x, z, systemRef.model,T);
    for iInt = 1:nInt
        integrand(iInt)  = integrandHandle{iInt}.intF(x,z,tau,t,systemRef);
    end
end