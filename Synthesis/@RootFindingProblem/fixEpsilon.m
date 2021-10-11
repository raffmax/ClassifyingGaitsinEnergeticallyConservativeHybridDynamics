function fixEpsilon(obj,epsLogical)
%fixEpsilon fixes/unfixes epsilon from the RootFindingProblem
%   epsLogical is a true or false and removes or adds epsilon to the
%   decision variables, respectively.
    obj.fixedEpsilon = epsLogical;
    obj.updateNumDecVar();
end
