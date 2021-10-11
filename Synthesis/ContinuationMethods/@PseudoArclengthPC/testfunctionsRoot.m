function testfunctionsRoot(obj, LCBefore, LCAfter)
    %testfunctionsRoot tests whether the current LC is a bifurcation
    % or if it lies after a bifurcation
    % Input: 
    %   LCBefore
    %   LCAfter 
    % Flags:
    %   1 - root of augmented jacobian
    %   2 - root of a time domain
  
    
    % get tanget spaces
    t1 = LCAfter.Sol.rfpData.jacobianAug(end,:)';
    t2 = LCBefore.Sol.rfpData.jacobianAug(end,:)';
    
    d = obj.Direction;
    
    % - 1 - Allgower: (8.1.17) Jumping Over A Bifurcation Point.
    if t1'*t2 < 0
        if obj.dispBifurcation
            disp(['There is a bifurcation between LC ',LCBefore.Name,' and LC ',LCAfter.Name,'!'])
        end
        % reverse orientation of curve
        obj.Direction = -obj.Direction;
%         LCAfter.setAfterBifurcation(1); 
%         LCBefore.setBeforeBifurcation(1);

        % update direction in LCAfter as well 
        % following LCs will be constructed with the updated
        % obj.Direction
        solAfter                   = LCAfter.Sol;
        solAfter.rfpData.direction = obj.Direction; % see few lines above
        LCAfter.updateLCSol(solAfter);

%        if d == 1
%            LCAfter.setAfterBifurcation(1); 
%            LCBefore.setBeforeBifurcation(1);
%        else
%            LCAfter.setBeforeBifurcation(1); 
%            LCBefore.setAfterBifurcation(1);
%        end

        LCAfter.setAfterBifurcation(1); 
        LCBefore.setBeforeBifurcation(1);
        
%         if direction == 1
%             LC1.setAfterBifurcation(1); 
%             LC2.setBeforeBifurcation(1);
%             
%             % update direction in LCAfter
%             % following LCs will be constructed with the updated
%             % obj.Direction
%             solAfter               = LC1.Sol;
%             solAfter.rfpData.direction = obj.Direction;
%             LC1.updateLCSol(solAfter);
%         else
%             LC1.setBeforeBifurcation(1);
%             LC2.setAfterBifurcation(1);
%             
%             % update direction in LCAfter
%             solAfter               = LC2.Sol;
%             solAfter.rfp.direction = obj.Direction;
%             LC2.updateLCSol(solAfter);
%         end
    end
    
    % - 2 -  check if a timeDomain is zero or changed its sign
    tdA_tdB = LCAfter.Sol.tDomain.*LCBefore.Sol.tDomain;
    if ~all(LCAfter.Sol.tDomain) % if not all entries are non-zero
        LC.setBifurcation(2); 
    elseif ~all(tdA_tdB >= 0) % if any timeDomain changed its sign
        rootIdx = find(tdA_tdB < 0); 
        LCAfter.setAfterBifurcation(2);
        LCBefore.setBeforeBifurcation(2);
        fprintf('Change in timeDomain(s) %d detected.', rootIdx)
    end
    

end