function nInt = getNumIntegrands(obj)
%getNumIntegrands counts all additional integrands
%   Returns the number of integrands defined in RootFindingProblem

nInt = 0;
if iscell(obj.functionals)
    for i = 1:length(obj.functionals)
       if isempty(obj.functionals{i})
       elseif obj.functionals{i}.int
           nInt = nInt+1;
       end
    end
end
end