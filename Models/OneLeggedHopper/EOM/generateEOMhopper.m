function model = generateEOMhopper(model)
%generateEOM – Creates executable handles of the model
%handles to mass matrix, coriolis, centrifugal,gravitational & external 
%forces, input mappings and collisions
%The EOM has the form:  M(q)u_dot = h(q, q_dot) + B(q)u + W(q, d_dot)lambda   
%Note: h(q, q_dot) includes tau_j, the passiv dynamic actuation (here via the spring)
%
% Syntax: [M_handle, h_handle, B_handle, W_handle, W_dot_handle, Wdyn_handle] = generateEOM(model)
%
% Inputs:       model  –   contains parameters of the model and the savepath
%                                for the optimized matlabfunctions
%
% Outputs:
% handles
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        University of Stuttgart, Institute for Nonlinear Mechanics
% Author:        Maximilian Raff
% email address: raff@inm.uni-stuttgart.de
% Website:       https://www.inm.uni-stuttgart.de/en
% Last revision: 12-July-2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Computing function handles for One Legged Hopper')
    %% --- Variables ---
    %  - generalized coor./vel. -
    syms x y l             real
    syms dx dy dalpha dl  real  
    alpha   = sym('alpha','real'); 
    
    epsilon = sym('epsilon','real');

    % - model parameters -
    m_t     = sym('m_t','real');    % mass of upper body / torso
    m_f     = sym('m_f','real');    % mass of foot / lower body
    g       = sym('g','real');      % gravity
    l_0     = sym('l_0','real');    % natural length of spring
    k_l     = sym('k_l','real');      % spring stiffness
    k_h     = sym('k_h','real');      % spring stiffness
    
    %  - discrete states -
    syms z_S z_F real  %Stance & Flight

    %% --- Vector Composition ---
    %  - generalized coor./vel. -
    q  = [x y alpha l]; %#ok<NODEF>
    dq = [dx dy dalpha dl];

    %  - discrete states -
    z = [z_S z_F];
    
    %  constraint functionals
    con_S = [x+l*sin(alpha);y-l*cos(alpha)];   % constraint during stance
    con_F = l-l_0; % constraint during flight
    
    %% --- Abbreviations ---
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       --- Lagrange ---                              %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% --- Kinematics ---
    % CoG-positions (from kinematics):
    Angle_f   = alpha;
    CoG_t     = [x;y];
    CoG_f     = CoG_t + l* [sin(Angle_f);-cos(Angle_f)];
    
    % CoG-velocities (computed via Jacobians):
    d_CoG_t   = jacobian(CoG_t,q)*dq.';
    d_CoG_f   = jacobian(CoG_f,q)*dq.';
    
    d_Angle_f   = jacobian(Angle_f,q)*dq.';
    
    %  - Potential Energy (due to gravity) -
    V_grav = CoG_t(2)*m_t*g + CoG_f(2)*m_f*g;
    %  - Potential Energy (due to springs) -
    V_spring  = 0.5*k_l*(l_0-l).^2 + 0.5*k_h*(0-alpha).^2;
    % Kinetic Energy:         
    T = 0.5 * (m_t     * sum(d_CoG_t.^2) + ...
               m_f     * sum(d_CoG_f.^2));
    
    V         = V_grav + V_spring;
              
    %% --- Lagrange 2 ---
    % only V_grav since the springs are dealt with via the contact matricies
    [M, h] = model.EulerLagrange(T, V_grav, q, dq); 
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    %  - Spring Forces -
    FspringL = k_l*(l_0-l); 
    FspringH = k_h*(0-alpha); 

    %% --- Projection Matricies ---

    B = [[0;0;1;0],...
        [ -(1-epsilon)*sin(alpha);(1-epsilon)*cos(alpha) ;0; epsilon ]];

    tau_J = B*[FspringH;FspringL];
                        
    [W_S, W_S_dot] = model.computeContactMatrix(con_S, q, dq);
    [W_F, W_F_dot] = model.computeContactMatrix(con_F, q, dq);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% --- Scaling with epsilon ---
    % load parameters for evaluation with subs()
    p = model.parameters;
    
    m          = p.m; %#ok<*NASGU>
    g          = p.g;
    l_0        = p.l_0;
    m_f        = p.m_f;
    m_t        = p.m_t;
    k_l        = p.k_l;
    k_h        = p.k_h;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- Trafo ---
    D      = diag([1, 1, 1/epsilon, 1/epsilon]);
    M      = simplify(expand(subs(D*M)));
    B      = [[0;0;1;0],[ -(1-epsilon)*sin(alpha)*z_S;(1-epsilon)*cos(alpha)*z_S;0 ; 1 ]];%simplify(expand(subs(D*B)));
    tau_J  = simplify(expand(subs(D*tau_J)));
    h      = simplify(expand(subs(D*h)));
    Wdyn_S = simplify(subs(D*W_S*epsilon));
    Wdyn_F = simplify(subs(D*W_F*epsilon));
    T      = subs(T);
    V      = subs(V);

    W_S     = simplify(subs(W_S));
    W_F     = simplify(subs(W_F));
    
    h  = h + tau_J;
                                          
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %% --- Compute derivatives ---
    %  transpose arguments
    q                 = q';
    dq                = dq';
    z                 = z';
    
    x                 = [q;dq];
    
    %  define n-dimensional vector
    vec_n       = sym('vec_n', [4, 1]);

    lambda1      = sym('lambda1', [1, 1]);
    lambda2      = sym('lambda2', [2, 1]);
    tau          = sym('u', [2 1]);
    
    symCell = {{epsilon, q, dq, x, z};...
               T;...
               V;...
               {M,vec_n};...
               h;...
               {B,tau};...
               {{Wdyn_S,lambda2},{Wdyn_F,lambda1}};... % Wdyn
               {{W_S,vec_n},{W_F,vec_n}};... % W
               {{W_S_dot,vec_n},{W_F_dot,vec_n}}... % Wdot
               };
           
    % generate matlab handels 
    model.generateMatlabHandels(symCell);
    
end