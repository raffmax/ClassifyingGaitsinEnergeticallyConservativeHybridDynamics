function options = setOptionsfsolve(obj)
%setOptionsfsolve sets options for fsolve
%   reads RootFindingProblem.options and sets fsolve options accordingly
options = optimoptions('fsolve');
fn = fieldnames(options);
for k=1:numel(fn)
    if isfield(obj.options,fn{k})
        options.(fn{k}) = obj.options.(fn{k});
    elseif strcmp(fn{k},'SpecifyObjectiveGradient')
        options.(fn{k}) = obj.options.Grad1;
    end
end
end

