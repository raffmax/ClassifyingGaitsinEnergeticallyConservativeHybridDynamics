function [exitflag,stepsPC] = runPC(obj,conPar,conPar_target)
    %runPC executes Predictor-Corrector Continuation
    % conPar = continuation parameter (lambda)
    % decVar = decision variable (x)
    % conPar_target = target of continuation
    if nargin <3
        conPar_target = [];
    end
    % check if epsilon gradients are disabled
    if ~obj.RFP.fixedEpsilon
        % this is necessary to seperate decVar and conPar
        obj.RFP.fixEpsilon(true);
    end

    % get decison variable without the continuation variable of current Limit Cycle
    % e.g. x = decVar and lambda = E0
    LCsol  = obj.LC.Sol;
    decVar = obj.RFP.getDecisionVariable(LCsol);

    % update epsilon in model reference
    obj.RFP.model.setEpsilon(LCsol.epsilon);
    % get full continuation variable (u in H(u))
    var = obj.getConVar(decVar,conPar);
    % get jacobian: fun = H, jacobian = HPrime(u) (Allgower2003)
    [fun,jacobian] = obj.getH(var);
    
    % compute tangent
    t = getTangent(jacobian);
    
    % get continuation orientation/direction
    if isempty(conPar_target)
        d = obj.Direction;    
    else
        if conPar < conPar_target
            d = 1;   
            if sign(t(obj.idxConPar))==1
                obj.Direction = 1; % orientation
            else
                obj.Direction = -1; % orientation
            end 
        else
            d = -1; 
            if sign(t(obj.idxConPar))==-1
                obj.Direction = 1; % orientation
            else
                obj.Direction = -1; % orientation
            end 
        end
    end
    
    % update augmented Jacobian in LC.Sol
    LCsol.rfpData.jacobianAug = [jacobian;t'];
    % update direction of initial LC
    LCsol.rfpData.direction = obj.Direction;
    obj.LC.updateLCSol(LCsol);
    % set (possibly) new controller
    obj.LC.setController(obj.RFP.controller)
    
    iter     = 0;
    stepsPC = zeros(obj.MaxIterations,1);

    if obj.dispIter
        disp('Iteration:');
    end
    
    
    while iter < obj.MaxIterations
        if ~isempty(conPar_target)
           if d*conPar+obj.StepTolerance >= d*conPar_target
               break
           end
        end
        if obj.dispIter
            if mod(iter, 50)==49 % print a '.' for each iter or '50', '100' etc. 
                fprintf('%d', iter+1); % iter starts with 0
            else
                fprintf('.');
            end
        end   
        
        % get full continuation variable (u in H(u))
        var = obj.getConVar(decVar,conPar);

        %% Predictor
        [varP,funP,jacobianP,tP,exitflagP,h,aimOnTarget] = obj.predictor(d,conPar_target,var,fun,t);
           
        %% Corrector
        [solC,exitflagC,funC,jacobianC,tC] = obj.corrector(varP,funP,jacobianP,tP,d,aimOnTarget,conPar_target);
        
        %% add Limit Cycle
        if exitflagC > 0
            % create and add new limit cycle
            newLC  = obj.addLC(solC,d);
            newLC.setLabel(obj.LC.Label);
            % update Limit Cycle and conPar, decVar, fun, t
            obj.LC = newLC;
            conPar = obj.getConPar(obj.LC.Sol);
            decVar = solC.rfpData.decVar;
            fun    = funC;
            t      = tC;
            % check for bifurcation
            if d == 1
                obj.testfunctionsRoot(newLC.Prev, newLC);
            else
                obj.testfunctionsRoot(newLC.Next, newLC);
            end
            obj.Tree = newLC.updateLabel(obj.Tree);
        else
            warning('The corrector-step was not successful!')
            break
        end
                
        % update
        iter = iter+1;
        stepsPC(iter) = h;
    end
    % resize PC_steps
    stepsPC = stepsPC(1:iter);
    if iter < obj.MaxIterations
        if iter == 0
            exitflag = 1;
        else   
            exitflag = exitflagC;
        end
    elseif isempty(conPar_target)
        exitflag = 1;
    else
        warning('Continuation stopped prematurally')
        exitflag = 0;
    end
    
    if obj.dispIter
        fprintf('\n') % end output of iterations
    end
end