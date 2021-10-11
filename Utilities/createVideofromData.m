% clear
% close all
% clc

% running motion
motion = [0 , 1 , 0 , 0 , 0 , 0;...
          0 , 0 , 0 , 0 , 1 , 0];

% get data 
load Motion_E5_ContDavid.mat
data = path_motion_trajectories; 

%% postprocessing
for i=1:length(data)
    xOffset = data{i}.states{1,3}(1,end);
    data{i}.states{1,4}(1,1:end) = data{i}.states{1,4}(1,1:end)+xOffset;  
    data{i}.states{1,5}(1,1:end) = data{i}.states{1,5}(1,1:end)+xOffset;  
    data{i}.states{1,6}(1,1:end) = data{i}.states{1,6}(1,1:end)+xOffset;  
end

%% Epsilon values
epsilon1 = linspace(0,1,length(data));
epsilon2 = epsilon1;

%% Generate Animations
%  - load in model -
%  - load in model -
full3D = true; % change the 3D lighting effects
borders = false; % show edges for torso and legs
COG = true; % show COG markers
updateEpsilon = true;

% Colours [R,G,B] for torso;legs;springs
colours = [230, 230, 230; 160, 160, 160; 255, 0, 0]/255; 

% time-steps
dt = 2e-2;

Biped_graphic = PlanarBiped_Conservative(full3D,COG, borders,colours,epsilon1(1), epsilon2(1),updateEpsilon);

% --- create animation for this data case ---
F_complete = [];
numDomains = size(motion,2);
offset = mod(length(data)-2,numDomains);
% get new epsilon range
epsRange = [0,linspace(epsilon1(1),epsilon1(end),length(epsilon1)-1),1];
for i=1:length(data)
    
    iSE = 0; %index of start/end
    z = [];
    t = [];
    x = [];
    if i==1 || i==length(data)
        % start/end show whole stride
        if i==1
            for domains_= 1:numDomains
                tData = data{i}.times{domains_};
                tdomain = tData(1):dt:tData(end-1);
                x = [x;interp1(tData,data{i}.states{domains_}',tdomain)];
                t = [t,tdomain];
                z = [z,repmat(motion(:, domains_),[1,length(tdomain)])];
            end
            T = repmat(2*data{i}.times{3}(end),1,length(t));
        else
           % complete stride 
           if domains_== 6
               domainRange = 1:numDomains;
           else
               domainRange = [domains_+1:numDomains,1:numDomains,1:numDomains];
           end
           for domains_= domainRange
                tData = data{i}.times{domains_};
                tdomain = tData(1):dt:tData(end-1);
                x = [x;interp1(tData,data{i}.states{domains_}',tdomain)];
                t = [t,tdomain];
                z = [z,repmat(motion(:, domains_),[1,length(tdomain)])];
           end
           T = repmat(2*data{i}.times{3}(end),1,length(t));
        end
    else
        domains_= numDomains-mod((length(data)-offset-i-1),numDomains);
        tData = data{i}.times{domains_};
        tdomain = tData(1):dt:tData(end-1);
        x = [x;interp1(tData,data{i}.states{domains_}',tdomain)];
        t = [t,tdomain];
        z = [z,repmat(motion(:, domains_),[1,length(tdomain)])];
        T = linspace(2*data{i-1}.times{3}(end),2*data{i}.times{3}(end),length(t));
    end
    F(length(t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    epsilon = linspace(epsRange(i),epsRange(i+1),length(t));
    for j= 1:length(t)
            state_ = x(j,:);
            z_     = z(:,j);
            % set actuation 
            u = [0;0;0;0];

            update(Biped_graphic, state_', u,z_, epsilon(j), epsilon(j));
            txtEps = text(state_(1)-1.5,0,2.1,['$\varepsilon= ',num2str(round(epsilon(j),2)),'$'],'FontSize',18,'Color','blue','Interpreter','latex');
            txtE   = text(state_(1)+.4,0,2.1,'$E_0=5.0~m_ol_o g$','FontSize',18,'Color','black','Interpreter','latex');
            txtT   = text(state_(1)+.4,0,1.9,['$T~= ',sprintf('%.2f',T(j)),'~\sqrt{l_o/g}$'],'FontSize',18,'Color','black','Interpreter','latex');
            % plot rectangle
            x1 = state_(1)-1.55;
            x2 = state_(1)-.92;
            x3 = x2;
            x4 = x1;
            y  = 0;
            z1 = 2.2;
            z2 = 2.2;
            z3 = 2;
            z4 = 2;
            hold on
            rectangle = plot3( [x1 x2 x3 x4 x1], [y y y y y], [z1 z2 z3 z4 z1],'b');
            % append
            F(frame_counter) = getframe(gcf);
            delete(txtEps)
            delete(txtT)
            delete(txtE)
            delete(rectangle)
            frame_counter = frame_counter + 1;
     end
     F_complete = [F_complete, F];
     clearvars F % s.t. F can be resized
     if iSE==length(data)
        % end animation
        break
     end
end

% -- make gif and video --
% video, note: more datapoints are needed for a good video
video_gif_ID = num2str(1);
    v = VideoWriter(strcat(pwd, '\AnimatioonBiped_', video_gif_ID, '.avi'));
    v.FrameRate = 40;
    v.Quality = 100;
    open(v)
    writeVideo(v,F_complete)
    close(v)
clearvars v 

clearvars F_complete

