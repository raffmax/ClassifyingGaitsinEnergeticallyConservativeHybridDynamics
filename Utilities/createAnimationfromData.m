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

% select data point
data = data{1};

%% Epsilon values
epsilon1 = 0;
epsilon2 = 0;

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
dt = 7e-2;

Biped_graphic = PlanarBiped_Conservative(full3D,COG, borders,colours,epsilon1, epsilon2,updateEpsilon);

% --- create animation for this data case ---
F_complete = [];
for domains_= 1:length(data.times)
    z = motion(:, domains_);
    tData = data.times{domains_};
    t = tData(1):dt:tData(end-1);
    x = interp1(tData,data.states{domains_}',t);
    F(length(t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    for i=1:length(t)
        state_ = x(i,:);
        % set actuation 
        u = [0;0;0;0];

        update(Biped_graphic, state_', u, z);
        % append
        F(frame_counter) = getframe(gcf);
        frame_counter = frame_counter + 1;
    end
    F_complete = [F_complete, F];

    clearvars F % s.t. F can be resized
end


% -- make gif 
% video, note: more datapoints are needed for a good video
video_gif_ID = num2str(1);

%  - create new folder to store the results -
result_folder_path = char(strcat(pwd, '\AppData\', 'Animation'));
if ~exist(result_folder_path, 'dir')
	mkdir(result_folder_path)
end

movie2gif(F_complete, strcat(result_folder_path, '\GIFBiped', video_gif_ID, '.gif'),'LoopCount', inf,'DelayTime', 0)

clearvars F_complete

