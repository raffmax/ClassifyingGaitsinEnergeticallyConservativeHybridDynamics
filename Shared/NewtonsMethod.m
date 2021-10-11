function[x,fval,exitflag,output,jacobian] = NewtonsMethod(fun,x,options,varargin)
%NewtonsMethod implements Newton's method.
%   NewtonsMethod can be called by the corrector-step of a continuation
%   method or an arbitrary root-finding-problem
%   see also getJacobianFD

if nargin>3
    f        = varargin{1};
    jacobian = varargin{2};
    t        = varargin{3}; % nullspace of jacobian / tangent space
    
    funCounter = 0; % counter for function evaluations
    
    if options.Grad1
        useFD = false; % use finite differences
    else
        useFD    =  true;
        stepSize = options.FiniteDifferenceStepSize;
        % get dimensions of jacobian
        nF    = length(f);
        nX    = length(x);
    end
    
    if ~isempty(options.idxConPar)
        % create zero vector and set the entry at the lambda position
        % equal to 1
        lambdaPos                       = zeros(1,length(x));
        lambdaPos(options.idxConPar)    = 1;
    else
        lambdaPos = [zeros(1,length(x)-1),1];
    end
    
else
    t = []; % tangent space is assumed to be empty
    
    % this option is only important in Continuation-Corrector
    options.aimOnTarget = false;
    
   if options.Grad1
       useFD = false;
       [f,jacobian] = fun(x);
       funCounter   = 1; % counter for function evaluations
   else
       useFD        = true;
       stepSize     = options.FiniteDifferenceStepSize;
       [f,jacobian] = getJacobianFD(fun,x,stepSize);
       % get dimensions of jacobian
       nF = length(f);
       nX = length(x);
       funCounter   = 1 + nX*nF; % counter for function evaluations
   end   
end

err  = max(abs(f));
iter = 0;
while err > options.FunctionTolerance && iter < options.MaxIterations  
   % newton step
   if options.aimOnTarget
        delta_x = -[jacobian;lambdaPos]\[f;0]; % fix lambda
   else
        delta_x = -[jacobian;t']\[f;zeros(size(t,2))]; % minimze distance to curve if size(t,2)=1 / equivalent to Moore-Penrose Inverse
   end

   if sum(abs(delta_x)) < options.StepTolerance 
       warning(['The Newton step is smaller than StepTolerance ',num2str(options.StepTolerance),'!'])
       break
   end

   % update
   x = delta_x + x;
   if useFD
       [f,jacobian] = getJacobianFD(fun,x,stepSize);
       funCounter   = funCounter + 1 + nX*nF; 
   else
       [f,jacobian] = fun(x);
       funCounter   = funCounter + 1;
   end
   
   % compute new tangent vector
   t = getTangent(jacobian);

   iter = iter + 1;
   err  = max(abs(f));       
end  
fval   = f;
if err < options.FunctionTolerance
   exitflag = 1;
else
   exitflag = -1;
end
output.t          = t;
output.iterations = iter;
output.funcCount  = funCounter;

end