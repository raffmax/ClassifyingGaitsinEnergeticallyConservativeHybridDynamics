classdef HopperParameters < handle
    % HopperParameters contains all the mechanical parameters for the 1D pogo stick.
    %   This class exists in order to not clutter the functionality in
    %   actual pogo stick class.
    %
    %   Input:  
    %      coefRes       coefficient of restitution
    %
    %   Methods:
    %      PogoStickParameters   The constructor for this class
    %
    %   Example:
    %      parameters = PogoStickParameters2D(coefRes);
    %
    %   See also PogoStick1D, MechanicalModel
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 12-July-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        coefRes            % coefficient of restitution                                                  
         
        % - pogo stick parameters -
        m       % total mass
        m_t     % mass of torso / upper body relative to m
        m_f     % mass of foot  / lower body relative to m
        theta_t % torso inertia 
        g       % gravity
        l_0     % natural length of spring
        k_l     % leg spring stiffness
        k_h     % hip spring stiffness    
    end
    
%     properties (Constant)
%         g                           = 1                                         % normalized gravity
%         l_0                         = 1                                         % normalized leglenght
%     end
    
    methods
        % - constructor - 
        function this = HopperParameters(coefRes,epsilon)
            % set mass and actuation scaling values
            this.coefRes = coefRes;
            if nargin <2
                epsilon = sym('epsilon','real');
            end
            
            % physical values
            this.m          = this.param_m(epsilon);  % total mass
            this.m_f        = this.param_m_f(epsilon);
            this.m_t        = this.param_m_t(epsilon);
            this.g          = this.param_g(epsilon);
            this.l_0        = this.param_l_0(epsilon);
            this.theta_t    = this.param_theta_t(epsilon);
            this.k_l        = this.param_k_l(epsilon);
            this.k_h        = this.param_k_h(epsilon);
        end
        function m = param_m(~,~)
            m = 1;
        end
        function m_t = param_m_t(this,~)
            m_t = this.m-this.m_f;
        end
        function m_f = param_m_f(this,epsilon)
            m_f = .05*this.m*epsilon;
        end
        function g = param_g(~,~)
            g = 1;
        end
        function l_0 = param_l_0(~,~)
            l_0 = 1;
        end
        function theta_t = param_theta_t(this,epsilon)
            l_0 = this.l_0;
            r   = l_0; % radius of gyration
            m_t = this.param_m_t(epsilon);
            theta_t = m_t*r^2;
        end
        function k_l = param_k_l(~,~)
            k_l = 40;
        end
        function k_h = param_k_h(this,epsilon)
            m_f = this.param_m_f(epsilon);
            g    = this.g;
            l_0  = this.l_0;
            k_h = 5*m_f*g*l_0;
        end
    end
end

