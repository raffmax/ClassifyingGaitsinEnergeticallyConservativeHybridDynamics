function [LCroot, flag] = approximateRoot(obj, LCBefore, LCAfter, type)
    %APPROXIMATEROOT find solutions that is an approximate a zero
    %point of a given functional
    %
    %   Based on Allgower chap. 9.2 -  Calculating Zero Points f(c(s))=0
    %
    %   Input:
    %       arcPC   - for using a predictor/corrector method
    %       LCBefore - LC before a zero point of a functional is detected (e.g.
    %                 after a bifurcation)
    %       LCAfter - LC after a zero point of a functional is detected (e.g.
    %                 after a bifurcation)
    %       type (string) - functional f(c(s)) which zero point should be
    %                 approximated. Implemented so far: simple bifurcation 
    %                 points/ 'detJacAug'
    %   Output:
    %       LCroot - Limit Cycle with solution of approximated zero point
    %       flag - 1 - loop terminated successfully with |h| <= hmin
    %              2 - loop terminated successfully with |f(c(s))| <= fmin
    %              -1 - iterations exceed maxIter
    
    %% setup
    
    obj.dispBifurcation = false; % turn off notification for bif temporarily
    
    direction        = obj.Direction; % save to reset after approximation
    obj.Direction    = LCBefore.Sol.rfpData.direction;
    LCAfterOriginal  = LCAfter; 
    LCBeforeOriginal = LCBefore; 
    solC             = LCBefore.Sol;
    
    % distance of Solution of LCBefore and LCAfter
    h = obj.getSolDistance(LCAfter.Sol,LCBefore.Sol)*1.5;
    %     h = abs(obj.getSolDistance(LCAfter.Sol,LCBefore.Sol))*obj.Direction; % kept here for debugging
    %     h = obj.getSolDistance(LCAfter.Sol, LCBefore.Sol); 
    
    % get root functional
    [rootFcn, typeRoot, opt]   = obj.getRootFindingFunctional(type);
    fvalBefore = rootFcn(LCBefore);
    fvalAfter  = rootFcn(LCAfter);
    
    if fvalAfter*fvalBefore > 0 % test 
        error("Something is wrong. Couldn't find the LC before the root of the functional %s", type);
    end
    fprintf('Approximating the root of %s between LC %s and  LC %s ...\n',...
        type, LCBefore.Name, LCAfter.Name);
  
    
    % stopping criteria
    h_min         = obj.RFP.options.StepTolerance; % stopping criteria Allgower
    fcn_min       = obj.RFP.options.FunctionTolerance; % tolerance for fval
    maxIter       = 50;
    if isfield(opt, 'stopCriteria') 
        fcnStop = opt.stopCriteria; % alternative stop criteria (for fcn_min), e.g smallest eigenvalue for detJacAug
    else
        fcnStop = @(fval) abs(fval);
    end
    
    F = zeros(maxIter,1); % Functional
    H = zeros(maxIter,1); % steplength

    %% loop until convergence, tolerance or maxIter
    iter = 1;
    while abs(h)>h_min && iter<=maxIter && min(fcnStop(fvalBefore),fcnStop(fvalAfter))>fcn_min
        F(iter) = fvalAfter;
        
        % set LCAfter as current LC in the continuation
        obj.LC = LCAfter;
        
        %% predictor with newton type steplength adaption for root
        [varP,funP,jacobianP,tP,~,h,aimOnTarget] = predictorRoot(obj, LCBefore, LCAfter, rootFcn, h);
        H(iter) = h;
        
        %% corrector
        [solC,exitflagC,~,~,~] = obj.corrector(varP,funP,jacobianP,tP,sign(h), aimOnTarget, []);
        
        %% create Limit Cycle
        % add only the final LCBif
        if exitflagC > 0
            % create new limit cycle
            newLC  = LimitCycle('to be deleted',solC,obj.RFP); % temporary limit cycle to further approximate the root
            obj.LC = newLC;
                            
            % add direction of prev. LCAfter to rfpData
            d = LCAfter.Sol.rfpData.direction;
            solC.rfpData.direction  = d;
            newLC.updateLCSol(solC);

            % set new LC After
            fvalBefore = fvalAfter;
            fvalAfter = rootFcn(newLC);
            LCBefore = LCAfter;
            LCAfter = newLC;
            
            if obj.dispIter
                if mod(iter, 50)==49 % print a '.' for each iter or '50', '100' etc.
                    fprintf('%d', iter+1); % iter starts with 0
                else
                    fprintf('.');
                end
            end
            
            iter = iter+1;
        end
        
    end
    
    if obj.dispIter
        fprintf('\n');
    end
    
    % output stopping criteria
    if abs(h)<=h_min
        flag = 1;
        fprintf('...terminated successfully with exit flag %i (convergence) \n', flag);
    elseif min(abs(fvalBefore),abs(fvalAfter))<=fcn_min
        flag = 2;
        fprintf('...terminated successfully with exit flag %i (function tolerance for root) \n', flag);
    else
        flag = -1;
        warning('... ended prematurly because maxIterations were reached.');
    end
    
    % set LC root if approximation was successful
    if flag > 0
%         LCstart = LCBeforeOriginal.getLCstart();
%         d = LCstart.Sol.rfpData.direction;
%         if d ~= obj.Direction
%             LCroot = obj.addLC(solC,d,LCBeforeOriginal);
%         else
%             LCroot = obj.addLC(solC,d,LCAfterOriginal);
%         end
%         LCroot = obj.addLC(solC,d,LCBeforeOriginal);
        if LCBeforeOriginal.Next == LCAfterOriginal
            LCroot = obj.addLC(solC,1,LCBeforeOriginal);
        else
            LCroot = obj.addLC(solC,-1,LCBeforeOriginal);
        end
        LCroot.setBifurcation(typeRoot);
        obj.Tree = LCroot.updateLabel(obj.Tree);
        obj.Bifurcations = [obj.Bifurcations LCroot];
        LCAfterOriginal.clearAfterBifurcation(1);
        LCBeforeOriginal.clearBeforeBifurcation(1);
    else
        LCroot = LimitCycle.empty();
    end
    
    obj.Direction = direction;
    
    % as sanitiy check; can be removed
    %     figure;
    %     subplot(1,2,1);
    %     plot(D, '.-');
    %     hold on;
    %     plot(D.*sign(H), '.-')
    %     xlabel('iterations');
    %     ylabel('fuctional of zero point');
    %
    %     subplot(1,2,2);
    %     stairs(H);
    %     xlabel('iterations');
    %     ylabel('steplength');
    
    obj.dispBifurcation = true; % turn on notification for bif again
    
end

