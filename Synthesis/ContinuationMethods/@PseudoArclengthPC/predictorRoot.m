function [uP,H_uP,Hprime_uP,tP,exitflagP,h,aimOnTarget] = predictorRoot(obj,LCBefore, LCAfter, rootFcn, h)
    % A euler predictor with a secant method for steplength adaptation
    % to approximate zero points on the solution curve
    % see Allgower[2003] chapter 9.2 - Calculating Zero Points
    % returns u_Predictor
    
    fcnBefore = rootFcn(LCBefore);
    fcnAfter = rootFcn(LCAfter);
    
    % predict wanted zero point at the zero crossing of the secant between
    % varAfter and varBefore 
    conParAfter = obj.getConPar(LCAfter.Sol);
    decVarAfter = LCAfter.Sol.rfpData.decVar;
    
    if obj.idxConPar > length(decVarAfter)
        u = [decVarAfter;conParAfter]; % conPar needs to be appended
    else
        u = decVarAfter; % conPar already included
    end
    
    
    % distance between solutions; direction points over root from LCAfter
    %     h0 = abs(obj.getSolDistance(LCAfter.Sol, LCBefore.Sol))*obj.Direction;
    
    % steplength to approximate root as zero point of secant
    % Allgower 9.2.3 (Newton Type Steplength Adaptation)
    % ratio h/h0 \in [-1,0]
    h = - fcnAfter/(fcnAfter-fcnBefore)*h;
    
    % tangent
    t_u = LCAfter.Sol.rfpData.jacobianAug(end,:)';
    
    % euler step
    v   = u + h*t_u;
    
    % result/output
    lambdaP        = v(obj.idxConPar);
    [H_v,Hprime_v] = obj.getH(v); % compute value and jacobian of the root function rootFunctionTDsingle();
    t_v            = getTangent(Hprime_v);
    
    uP        = v;
    H_uP      = H_v;
    Hprime_uP = Hprime_v;
    tP        = t_v;
    
    % values needed in the interface
    aimOnTarget = false;
    exitflagP = 1;
end