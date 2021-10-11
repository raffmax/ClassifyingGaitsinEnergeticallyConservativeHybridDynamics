function plotLagrangeMultipliers(lambdaMatrix,continuationVector)
figure('Name','Change of Lagrange Multipliers')
hold on
nL = size(lambdaMatrix,2);

for iL = 1:nL
    plot(continuationVector,lambdaMatrix(:,iL),'DisplayName',['$\lambda_{',num2str(iL),'}$'])
end
grid on
xlabel('Continuation Parameter','Interpreter','latex')
legend('Interpreter','latex')

end
