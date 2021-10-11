        function [LC1,LC2] = branchingOff(obj, LCBif)
        % this function finds the tangent directions in a simple
        % bifurcation point by approximating the bifurcation equation and
        % returns two limit cycles that lie in the so far unknown direction
        
            epsilon = obj.RFP.options.FiniteDifferenceStepSize*100; % this was just a guess, there is probably a better way to choose it
%             epsilon = 1e-6;
            
            conPar  = obj.getConPar(LCBif.Sol);
            decVar  = LCBif.Sol.rfpData.decVar;
            
            u_bar = obj.getConVar(decVar,conPar);           
            
            % see Allgower[2003] chapter 8.3 - Branching Off via
            % Bifurcation Equation
            HPrime  = LCBif.Sol.rfpData.jacobianAug(1:end-1,:);
            tau1    = LCBif.Sol.rfpData.jacobianAug(end,:)';          
            [tau2,e]= approxKernels(HPrime,tau1);
            hessian = approxBifEq(obj,epsilon,u_bar,tau1,tau2,e);
            [t1,t2] = approxTangents(hessian,tau1,tau2);
            
            
            % decide which tangent direction is the so far unknown
            % direction
            if ~isempty(LCBif.Next)
                t    = LCBif.Next.Sol.rfpData.jacobianAug(end,:)';
                if abs(t'*t1) < abs(t'*t2)
                    t = t1;
                else
                    t = t2;
                end
            end
     
            %% Predictor
            conPar_target = [];
            [H_u,~] = obj.getH(u_bar,false);
            d = 1;
            [varP1,funP1,jacobianP1,tP1,exitflagP1,h1,aimOnTarget] = obj.predictor(d,conPar_target,u_bar,H_u,t);
            [varP2,funP2,jacobianP2,tP2,exitflagP2,h2,aimOnTarget] = obj.predictor(d,conPar_target,u_bar,H_u,-t);

            %% Corrector
            [solC1,exitflagC1,funC1,jacobianC1,tC1] = obj.corrector(varP1,funP1,jacobianP1,tP1,d,aimOnTarget,conPar_target);
            [solC2,exitflagC2,funC2,jacobianC2,tC2] = obj.corrector(varP2,funP2,jacobianP2,tP2,d,aimOnTarget,conPar_target);

            [conPar1, name1] = obj.getConPar(solC1);
            [conPar2, name2] = obj.getConPar(solC2);
            
            % limitcycles on the new branch
            LC1    = LimitCycle([name1,'=',num2str(conPar1)],solC1,obj.RFP);
            LC2    = LimitCycle([name2,'=',num2str(conPar2)],solC2,obj.RFP);
            
            
            % add more rfp data to solC for locateBifurcation
            if obj.Direction*t'*tC1 > 0
                solC1.rfpData.direction = d;
                solC2.rfpData.direction = d;
                LCBif.addNext(LC1,LC2);
            else
                solC1.rfpData.direction = -d;
                solC2.rfpData.direction = -d;
                LCBif.addPrev(LC1,LC2);
            end
            
            LC1.updateLCSol(solC1);
            LC2.updateLCSol(solC2);

            
            obj.Tree = LC1.updateLabel(obj.Tree,LC2);
            
        end
        
        function [tau2,e] = approxKernels(Hprime,tau1)
             A          = [Hprime; tau1'];
             [tau2,~]   = eigs(A'*A, 1, 'SM');
             [e,~]      = eigs(Hprime*Hprime', 1, 'SM');
        end
        
        function hessian = approxBifEq(obj, epsilon, u_bar,tau1,tau2,e)
            %% would be better to have a function that only computes the value of f and 
            %not also the jacobian like getH(...) does
            
            % function evaluations for finite differences
            ge0     = e'*obj.getH(u_bar+epsilon*tau1,false);
            g00     = e'*obj.getH(u_bar,false);
            g_e0    = e'*obj.getH(u_bar-epsilon*tau1,false);
            g0e     = e'*obj.getH(u_bar+epsilon*tau2,false);
            g0_e    = e'*obj.getH(u_bar-epsilon*tau2,false);
            gee     = e'*obj.getH(u_bar+epsilon*tau1+epsilon*tau2,false);
            g_e_e   = e'*obj.getH(u_bar-epsilon*tau1-epsilon*tau2,false);
            ge_e    = e'*obj.getH(u_bar+epsilon*tau1-epsilon*tau2,false);
            g_ee    = e'*obj.getH(u_bar-epsilon*tau1+epsilon*tau2,false);
            
            % finite differences
            a11 = 1/epsilon^2*(ge0-2*g00+g_e0);
            a22 = 1/epsilon^2*(g0e-2*g00+g0_e);
            a12 = 1/(4*epsilon^2)*(gee+g_e_e-ge_e-g_ee);
            a21 = a12;
            
            hessian = [a11 a12; a21 a22];
        end
        
        function [t1,t2] = approxTangents(hessian,tau1,tau2)
        
            if hessian(1,1) ~= 0
                r = roots([hessian(1,1), 2*hessian(1,2), hessian(2,2)]);
                t1 = r(1)*tau1 + tau2;
                t2 = r(2)*tau1 + tau2;
            else
                r = roots([hessian(2,2), 2*hessian(1,2), hessian(1,1)]);
                t1 = tau1 + r(1)*tau2;
                t2 = tau1 + r(2)*tau2;
            end
            
            t1 = t1/norm(t1);
            t2 = t2/norm(t2);
        end
        