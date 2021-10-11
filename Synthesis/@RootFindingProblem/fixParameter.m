function fixParameter(obj,paramLogical)
%fixParameter fixes/unfixes parameter from the RootFindingProblem
%   paramLogical is a true or false and removes or adds parameter to the
%   evaluation of user-provide functionals, respectively.
    obj.fixedParameter = paramLogical;
end
