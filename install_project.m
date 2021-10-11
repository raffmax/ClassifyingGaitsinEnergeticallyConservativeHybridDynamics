%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   --- Project Installer ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Describtion:	Adds the paths for the ModelHomotopy libraries, scripts, and
% tests files.
%
% Output:       res     -     returns 1 if task was successful
%
% Example:      - ~
%
%
% Other m-files required:   none
% MAT-files required:       none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        University of Stuttgart, Institute for Nonlinear Mechanics
% Author:        Maximilian Raff
% email address: raff@inm.uni-stuttgart.de
% Website:       https://www.inm.uni-stuttgart.de/en
% Last revision: 23-March-2021

function res = install_project()
%% --- clean workspace ---
clear
close all
clc

%% --- add necessary directiroes to the pathdef ---
dir_ModelHomotopy = pwd;                                                    % get a handle on the current path
path(pathdef);                                                              % reset path to its default value
addpath(genpath(dir_ModelHomotopy));                                        % add necessary folders of the model on top of the path

%% delete AppData or remove it from matlab path
res = deleteAppData();
if res<0
    rmpath(genpath([dir_ModelHomotopy, filesep, 'AppData']));
end
%% --- return ---
res = 1;                                                                    % install successful
end