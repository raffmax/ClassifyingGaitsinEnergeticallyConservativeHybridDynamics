function removeFunctional(obj,idx)
%removeFunctional removes user-defined functionals for RootFindingProblem
%   idx  - position/index of functional
    if idx==1 && length(obj.functionals{:})==1
        obj.functionals = [];
        obj.fixedParameter(false);
    else
        obj.functionals{idx} = [];
    end
    % update
    obj.nInt    = obj.getNumIntegrands();
    obj.nCon    = obj.getNumConstraints();
    obj.updateNumDecVar();
end
