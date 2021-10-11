function plotTrajectory(trajData,stateID,newFigure,type,color)
%plotTrajectory plots chosen State/Actuation/Force over Time

if nargin <3
    newFigure = true;
end

if nargin <4
    type = 'x';
end

if nargin <5
    color = [0 0 1]; % blue
end

if newFigure
    figure
end
hold on
grid on
switch length(stateID)
    case 1
        for iDomain=1:length(trajData)
            plot(trajData(iDomain).t,trajData(iDomain).(type)(:,stateID),'color',color)
            scatter(trajData(iDomain).t(1),trajData(iDomain).(type)(1,stateID),'k','o','filled')
            scatter(trajData(iDomain).t(end),trajData(iDomain).(type)(end,stateID),'k','o')
            xlabel('$t$','Interpreter','latex')
        end
    case 2
        for iDomain=1:length(trajData)
            plot3(trajData(iDomain).t,trajData(iDomain).(type)(:,stateID(1)),...
                                      trajData(iDomain).(type)(:,stateID(2)),'color',color)
            xlabel('$t$','Interpreter','latex')
        end
end
end

