function plotEnergyTimeTrajectory(LC,trajData,z)
%plotEnergyTimeTrajectory plots system's energy over Time

nDomains = size(trajData,2);
LC.model.setEpsilon(LC.Sol.epsilon);
figure('Name','Energy-Trajectory')
hold on
for iDomain = 1:nDomains
    E = zeros(size(trajData(iDomain).t));
    for iE = 1:length(E)
        E(iE) = LC.model.SystemEnergy(trajData(iDomain).x(iE,:)',z(:,iDomain));
    end
    plot(trajData(iDomain).t,E)
end
xlabel('$t$','Interpreter','latex')
ylabel('$E(x)$','Interpreter','latex')
title(['$\epsilon =',num2str(LC.model.epsilon),'$'],'Interpreter','latex')
grid on
end