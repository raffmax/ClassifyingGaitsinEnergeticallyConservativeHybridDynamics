function setFunctional(obj,idx,type,varargin)
%setFunctional user-defined functionals for RootFindingProblem
%   idx  - position/index of functional
%   type - 'start','end' or pre-implemted functionals: e.g. 'E0', 'E_avg' 

    p = inputParser;
   
    % default: Energy constraint on E0
    z0            = obj.sequence.FlowMapOrder(:,1);
    default_F     = @(x,tau,c,rfp) rfp.model.SystemEnergy(x,z0)-c;
    default_F_x   = @(x,tau,c,rfp) rfp.model.Grad1.E_x(rfp.model.getEpsilon, x, z0);
    default_F_tau = @(x,tau,c,rfp) zeros(1,size(rfp.controller.getB,2));
    default_F_eps = @(x,tau,c,rfp) rfp.model.Grad1.E_eps(rfp.model.getEpsilon, x, z0);
    
    % cfg 
    addParameter(p, 'c',        'E0');          % string of constant parameter
    addParameter(p, 'cost',     false);         % is functional also used as objective/cost function?
    addParameter(p, 'int',      false);         % is functional obtained by user-defined integration?
    addParameter(p, 'F',        default_F);     % function value F at t=0 (start) or t=T (integral,end)
    addParameter(p, 'F_x',      default_F_x);   % its derivative w.r.t. x
    addParameter(p, 'F_tau',    default_F_tau); % its derivative w.r.t. tau
    addParameter(p, 'F_eps',    default_F_eps); % its derivative w.r.t. epsilon
    addParameter(p, 'F_T',      []);            % its derivative w.r.t. period T
    addParameter(p, 'F_int',    []);            % its derivative w.r.t. user-defined integral
    addParameter(p, 'intF',     []);            % function value F(x,tau) (e.g. running cost)
    addParameter(p, 'intF_x',   []);            % its derivative w.r.t. x
    addParameter(p, 'intF_tau', []);            % its derivative w.r.t. tau
    addParameter(p, 'intF_eps', []);            % its derivative w.r.t. epsilon

    % load pre-implemted functionals
    switch type
        case 'E0' % default
            newType = 'start';
            obj.fixedParameter(true);
        case 'E_avg'
            varargin = ...
            {'int'     , true,...
            'c'       , 'E_avg',...
            'intF'    ,@(x,z,tau,t,obj) obj.model.SystemEnergy(x,z),...
            'intF_x'  ,@(x,z,tau,t,obj) obj.model.Grad1.E_x(obj.model.getEpsilon,x,z),...
            'intF_tau',@(x,z,tau,t,obj) zeros(1,obj.model.nTau),...
            'intF_eps',@(x,z,tau,t,obj) obj.model.Grad1.E_eps(obj.model.getEpsilon,x,z),...
            'F'       ,@(x,int,tau,T,c,obj) int-c*T,...
            'F_x'     ,@(x,int,tau,T,c,obj) zeros(1,2*obj.model.nQ),...
            'F_tau'   ,@(x,int,tau,T,c,obj) zeros(1,obj.model.nTau),...
            'F_eps'   ,@(x,int,tau,T,c,obj) 0,...
            'F_int'   ,@(x,int,tau,T,c,obj) 1,...
            'F_T'     ,@(x,int,tau,T,c,obj) -c};
            newType = 'end';
            obj.fixedParameter(true);
        case 'xDot_avg'
            varargin = ...
            {'c'    , 'xDot_avg',...
             'F'    ,@(x,int,tau,T,c,obj) x(1)-c*T,...
             'F_x'  ,@(x,int,tau,T,c,obj) [1 zeros(1,2*obj.model.nQ-1)],...
             'F_tau',@(x,int,tau,T,c,obj) zeros(1,obj.model.nTau),...
             'F_eps',@(x,int,tau,T,c,obj) 0,...
             'F_T'  ,@(x,int,tau,T,c,obj) -c};
             obj.fixedParameter(true);
            newType = 'end';
        case 'tau^2'
            varargin = ...
            {'int'     , true,...
            'cost'    , true,...
            'intF'    ,@(x,z,tau,t,obj) tau'*tau,...
            'intF_x'  ,@(x,z,tau,t,obj) zeros(1,2*obj.model.nQ),...
            'intF_tau',@(x,z,tau,t,obj) 2*tau',...
            'intF_eps',@(x,z,tau,t,obj) 0,...
            'F'       ,@(x,int,tau,T,c,obj) int,...
            'F_x'     ,@(x,int,tau,T,c,obj) zeros(1,2*obj.model.nQ),...
            'F_tau'   ,@(x,int,tau,T,c,obj) zeros(1,obj.model.nTau),...
            'F_eps'   ,@(x,int,tau,T,c,obj) 0,...
            'F_int'   ,@(x,int,tau,T,c,obj) 1,...
            'F_T'     ,@(x,int,tau,T,c,obj) 0};
            newType = 'end';
        otherwise
            newType = [];
    end
    if ~isempty(newType)
        type = newType;
    end
    
    parse(p, varargin{:}); 

    if length(obj.functionals)>=idx
        if ~isempty(obj.functionals{:})
            if ~isempty(obj.functionals{idx}.type)
                warning(['overwrite functionals{',num2str(idx),'}! Remove functional first the next time!'])
                % delete functional
                obj.removeFunctional(idx);
            end
        end
    end
    obj.functionals{idx}.type = type;
    obj.functionals{idx}.c    = p.Results.c;
    obj.functionals{idx}.cost = p.Results.cost;
    obj.functionals{idx}.int  = p.Results.int;
    
    
    obj.functionals{idx}.F     = p.Results.F;
    obj.functionals{idx}.F_x   = p.Results.F_x; 
    obj.functionals{idx}.F_tau = p.Results.F_tau;
    obj.functionals{idx}.F_eps = p.Results.F_eps; 
    switch type
        case 'end' % handles of type @(x,int,tau,T,c,obj)
            % for int* the handle is @(x,z,tau,T,obj)
            obj.functionals{idx}.intF     = p.Results.intF;
            obj.functionals{idx}.intF_x   = p.Results.intF_x;
            obj.functionals{idx}.intF_tau = p.Results.intF_tau;
            obj.functionals{idx}.intF_eps = p.Results.intF_eps;
           
            obj.functionals{idx}.F_T   = p.Results.F_T; 
            obj.functionals{idx}.F_int = p.Results.F_int; 
        case 'start' % handles of type @(x,tau,c,obj)
        otherwise
            error(['The type ,',type,' is not valid!'])
    end
    
    % update
    obj.nInt    = obj.getNumIntegrands();
    obj.nCon    = obj.getNumConstraints();
    obj.updateNumDecVar();
end
