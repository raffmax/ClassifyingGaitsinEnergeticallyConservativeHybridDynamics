function F_complete = getAnimationHopper(data)
%getAnimationHopper Animates the periodic solution
%   adapted from https://bitbucket.org/ganzhenyu/bipedalslipmodel_swingleg

% Physics:
systParam.g         = 1;     % [*] gravity
% Parameter of the model
systParam.l_0       = 1;     % [*] uncompressed leg length
systParam.m_0       = 1;     % [m_0] total mass
systParam.k         = 20;    % [m_0*g/l_0] linear spring stiffness in the leg
systParam.kh        = 5;     % [m_0*g/l_0] linear spring stiffness in the leg
systParam.st        = 0;     % 
systParam.Llo       = 0;     % 

hopper_graphic = SLIP_Model_Graphics_AdvancedPointFeet(systParam);

F_complete = [];
for iDomain = 1:size(data,2)
    F(length(data(iDomain).t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    for i=1:length(data(iDomain).t)
        state_ = [data(iDomain).x(i,1:2),0,data(iDomain).x(i,3:6),0,data(iDomain).x(i,7:8)]';
        z      = data(iDomain).z(i,:)';
        
        update(hopper_graphic, state_');
        % append
        F(frame_counter) = getframe(gcf);
        frame_counter = frame_counter + 1;
    end
    F_complete = [F_complete, F];

    clearvars F % s.t. F can be resized
end
end

