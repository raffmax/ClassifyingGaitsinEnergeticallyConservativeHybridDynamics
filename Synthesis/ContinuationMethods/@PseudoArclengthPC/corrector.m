function [solC,exitflagC,funC,jacobianC,tC] = corrector(obj,varP,funP,jacobianP,tP,d,aimOnTarget,conPar_target)
    % see also RFP.NewtonMethod()
    
    opts = obj.RFP.options;
    opts.aimOnTarget = aimOnTarget;
    opts.idxConPar = obj.idxConPar;
    [varC,funC,exitflagC,outputC,jacobianC] = NewtonsMethod(@(x)obj.getH(x),varP,opts,funP,jacobianP,tP);
  
    conPar = varC(obj.idxConPar);
    
    if ~isempty(conPar_target) && ~aimOnTarget && d*conPar+obj.StepTolerance >= d*conPar_target
        % corrector went over the conPar_target value
        % fix conPar to conPar_target
        varP(obj.idxConPar)        = conPar_target;
        opts.aimOnTarget = true;
        [varC,funC,exitflagC,outputC,jacobianC] = NewtonsMethod(@(x)obj.getH(x),varP,opts,funP,jacobianP,tP);
    end
    
    % get solution struct for corrector
    solC = obj.getSolStruct(varC);

    tC     = outputC.t;
   
    if obj.idxConPar == length(varC)
        decVar = varC(1:end-1);
    else
        decVar = varC;
    end
    
    % add more data to solution struct
    solC   = obj.RFP.getOutputSolvers(solC,decVar,outputC,jacobianC,exitflagC);
    % add more rfp data to solC for locateBifurcation
    solC.rfpData.direction = obj.Direction;
end