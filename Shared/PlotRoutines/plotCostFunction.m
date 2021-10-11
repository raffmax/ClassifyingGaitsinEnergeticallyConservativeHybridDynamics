function plotCostFunction(cost,continuationVector)
figure('Name','Cost of Limit Cycle')

plot(continuationVector,cost)
grid on
xlabel('Continuation Parameter','Interpreter','latex')
ylabel('Cost','Interpreter','latex')

end