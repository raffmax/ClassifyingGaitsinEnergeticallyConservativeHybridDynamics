function plotStateTimeTrajectories(trajData,stateLabel)
%plotStateTimeTrajectories plots all states over Time
% is dependent on plotTrajectory

if nargin <2
    stateLabel = [];
end

nStates = size(trajData(1).x,2);

for i = 1:nStates
    subplot(2,nStates/2,i)
    plotTrajectory(trajData,i,false)
    if ~isempty(stateLabel)
        ylabel(stateLabel{i},'Interpreter','latex')
    end
end
end

