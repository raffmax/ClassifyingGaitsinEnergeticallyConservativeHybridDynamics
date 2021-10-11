function [xE,intE] = getGradientsFromSimData(obj,simData)
%getGradientsFromSimData Computes gradients
%    Compute gradients of end states xE and intE of flowMaps
%    w.r.t. initial state x0, domain times t, parameters p 

if obj.options.Grad1
    fMD     = simData.flowMapData;
    jMD     = simData.jumpMapData;
    nDomain = obj.nTime;
    
    timeBased = obj.controller.timeBased;

    % allocate space
    xE.x0  = fMD.xGrad1_E;  
    xE.xi  = fMD.xiGrad1_E; 
    xE.eps = fMD.epsGrad1_E;
    
    intE.x0  = fMD.int.xGrad1_E; 
    intE.xi  = fMD.int.xiGrad1_E; 
    intE.eps = fMD.int.epsGrad1_E; 

    for j = 2:nDomain % loop over domains
        i = j-1;

        % D_ji * PhiParam_i
        xPlus_xi       = jMD.xGrad1(:,:,i)*xE.xi(:,:,i);
        if ~obj.fixedEpsilon
            % D_ji * PhiParam_i + Dparam_ji
            xPlus_eps  = jMD.xGrad1(:,:,i)*xE.eps(:,:,i)...
                         + jMD.epsGrad1(:,:,i);
        end
        
        % PhiParam_j + Phi_j *( D_ji * PhiParam_i + Dparam_ji )
        xE.xi(:,:,j) = xE.xi(:,:,j) + xE.x0(:,:,j)*xPlus_xi;
        if ~obj.fixedEpsilon
            xE.eps(:,:,j) = xE.eps(:,:,j)+ xE.x0(:,:,j)*xPlus_eps;
        end
        
        % Phi_j * D_ji * Phi_i
        xE.x0(:,:,j) = xE.x0(:,:,j)*jMD.xGrad1(:,:,i)*xE.x0(:,:,i);
        
        if obj.nInt >0
            % similar to steps before
            % Note: we take partial derivatives over a sum of integrals
            intE.xi(:,:,j) = intE.xi(:,:,j) + intE.x0(:,:,j)*xPlus_xi...
                             + intE.xi(:,:,i);
            if ~obj.fixedEpsilon
                intE.eps(:,:,j) = intE.eps(:,:,j)+ intE.x0(:,:,j)*xPlus_eps...
                                  + intE.eps(:,:,i);
            end 
            intE.x0(:,:,j) = intE.x0(:,:,j)*jMD.xGrad1(:,:,i)*xE.x0(:,:,i)...
                             + intE.x0(:,:,i);
        end
    end

    % - end states gradients w.r.t. time intervals/tDomains, i.e. f_j -
    xE.t            = zeros(obj.nState,nDomain,nDomain);
    intE.t          = zeros(obj.nInt,nDomain,nDomain);
    Phi             = fMD.xGrad1_E;
    D               = jMD.xGrad1;
    if timeBased % non-autonomous
        for ix = 1:nDomain
            for it = 1:nDomain
                xE.t(:,it,ix)  = fMD.TGrad1_E(:,ix);
                if obj.nInt >0
                    intE.t(:,it,ix) = fMD.int.TGrad1_E(:,ix);
                    if ix>1
                        % add previous partial derivatives, since every
                        % integral of a domain is summed up.
                        intE.t(:,it,ix) = intE.t(:,it,ix)+intE.t(:,it,ix-1);
                    end
                end
                if ix > it
                    xE.t(:,it,ix) = xE.t(:,it,ix) + fMD.tGrad1_E(:,ix);
                    
                    % similar for user provided integrals
                    if obj.nInt >0
                        intE.t(:,it,ix) = intE.t(:,it,ix) + fMD.int.tGrad1_E(:,ix);
                    end
                elseif ix == it
                    xE.t(:,it,ix) = xE.t(:,it,ix) + fMD.f_E(:,ix);
                    % similar for user provided integrals
                    if obj.nInt >0
                        intE.t(:,it,ix) = intE.t(:,it,ix)+fMD.integrand_E(:,ix);            
                    end
                end
                if ix > 1
                    xE.t(:,it,ix) = xE.t(:,it,ix) ... 
                                    + Phi(:,:,ix)*D(:,:,ix-1)*xE.t(:,it,ix-1);
                    % similar for user provided integrals
                    if obj.nInt >0
                        intE.t(:,it,ix) = intE.t(:,it,ix) ...
                                          + fMD.int.xGrad1_E(:,:,ix)*D(:,:,ix-1)*xE.t(:,it,ix-1); 
                    end
                end                      
            end
        end
    else % autonomous
        xE.t(:, 1, 1)   = fMD.f_E(:,1);
        intE.t(:, 1, 1) = fMD.integrand_E(:,1);
        for ix = 2:nDomain
            for it = 1:ix % only proceed to the diagonal since future times do NOT affect past states
                if it < ix
                    xE.t(:, it, ix)   = Phi(:,:,ix)*D(:,:, ix-1)*xE.t(:, it, ix-1);

                    if obj.nInt >0
                        intE.t(:, it, ix) = fMD.int.xGrad1_E(:,:,ix)*D(:,:, ix-1)*xE.t(:, it, ix-1)...
                                                + intE.t(:, it, ix-1);
                    end
                else
                    xE.t(:, it, ix)   = fMD.f_E(:,ix); % f_end of domain
                    if obj.nInt >0
                        intE.t(:, it, ix) = fMD.integrand_E(:,ix)+intE.t(:, it, ix-1);
                    end
                end
            end
        end
    end
else
    xE   = [];
    intE = [];
end
end

