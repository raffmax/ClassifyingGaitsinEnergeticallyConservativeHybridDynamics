function setSolverOptions(obj,options)
%setSolverOptions updates solver options

% set solver options
default_options.FunctionTolerance         = 1e-7;
default_options.MaxIterations             = 1000;
default_options.FiniteDifferenceStepSize  = 1e-8;
default_options.StepTolerance             = 1e-10;
default_options.MaxFunctionEvaluations    = 100*obj.nDecVar;
default_options.Grad1                     = false;
default_options.Grad2                     = false;
fn = fieldnames(default_options);
for k=1:numel(fn)
    if isfield(options,fn{k})
        % set new value for field provided by options
        obj.options.(fn{k}) = options.(fn{k});
    elseif ~isfield(obj.options,(fn{k}))
        % if field doesn't exist yet, use default
        obj.options.(fn{k}) = default_options.(fn{k});
    end
end

end

