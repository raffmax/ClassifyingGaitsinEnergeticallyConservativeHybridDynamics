function generalSetup(obj, varargin)
    % default values for the configuration
    default_fixedStates      = [];
    default_periodicStates   = (1:2*obj.model.nQ);
    default_periodicStatesHomotopy = [];
    default_fixedParameter   = true;
    default_fixedEpsilon     = true;
    default_optimalControl   = false;
    default_useLagrangeMultipliers = true;
    default_flipOperator     = eye(2*obj.model.nQ);
    default_controller       = ConservativeController(obj.model, 0);
    default_simulation       = Simulation();

    default_integrationType  = 'TD';
    default_shootingType     = 'single';
    default_options          = []; % will be set after parser


    p = inputParser;
    % cfg 
    addParameter(p, 'fixedStates',     default_fixedStates);
    addParameter(p, 'periodicStates',  default_periodicStates);
    addParameter(p, 'periodicStatesHomotopy',  default_periodicStatesHomotopy);
    addParameter(p, 'flipOperator',    default_flipOperator);
    addParameter(p, 'fixedParameter',  default_fixedParameter);
    addParameter(p, 'fixedEpsilon',    default_fixedEpsilon);
    addParameter(p, 'optimalControl',  default_optimalControl);
    addParameter(p, 'useLagrangeMultipliers',  default_useLagrangeMultipliers);
    addParameter(p, 'controller',      default_controller);
    addParameter(p, 'simulation',      default_simulation);

    addParameter(p, 'integrationType', default_integrationType);
    addParameter(p, 'shootingType',    default_shootingType);
    addParameter(p, 'options',         default_options);

    parse(p, varargin{:}); 

    obj.fixedStates        = p.Results.fixedStates;
    obj.freeStates         = default_periodicStates(~ismember(default_periodicStates,p.Results.fixedStates));
    obj.periodicStates     = p.Results.periodicStates;
    obj.periodicStatesHomotopy = p.Results.periodicStatesHomotopy;
    obj.flipOperator       = p.Results.flipOperator;
    obj.fixedParameter     = p.Results.fixedParameter;
    obj.fixedEpsilon       = p.Results.fixedEpsilon;
    obj.optimalControl     = p.Results.optimalControl;
    obj.useLagrangeMultipliers = p.Results.useLagrangeMultipliers;
    obj.controller         = p.Results.controller;
    obj.simulation         = p.Results.simulation;
    obj.simulation.setSystemRef(obj); % add reference
    obj.integrationType    = p.Results.integrationType;
    obj.shootingType       = p.Results.shootingType;
    
    obj.nState     = 2*obj.model.nQ;
    obj.nTime      = obj.getNumDomainTime();
    obj.nCtrl      = obj.controller.nXi;
    obj.nFreeState = length(obj.freeStates);
    obj.nCon       = obj.getNumConstraints();
    obj.nInt       = obj.getNumIntegrands();
    obj.nDecVar    = obj.getNumDecVar();
    
    % set solver options
    obj.setSolverOptions(p.Results.options)
    
    % set default structure for functionals
    if obj.fixedParameter
        obj.setFunctional(1,'E0')
    end
end