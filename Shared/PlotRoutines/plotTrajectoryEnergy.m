function plotTrajectoryEnergy(trajData,model)
%plotTrajectory plots chosen State/Actuation/Force over Time
figure
hold on
    for iDomain=1:length(trajData)
        E = zeros(size(trajData(iDomain).t));
        for i = 1:length(E)
            E(i) = model.SystemEnergy(trajData(iDomain).x(i,:)', trajData(iDomain).z(i,:)');
        end
        plot(trajData(iDomain).t,E,'b')
        scatter(trajData(iDomain).t(1),E(1),'k','o','filled')
        scatter(trajData(iDomain).t(end),E(end),'k','o')
        xlabel('$t$','Interpreter','latex')
        ylabel('$E$','Interpreter','latex')
    end
end

