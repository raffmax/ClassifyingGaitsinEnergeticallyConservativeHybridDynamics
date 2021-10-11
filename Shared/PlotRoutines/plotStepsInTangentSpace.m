function plotStepsInTangentSpace(continuationVector,PC_steps)
figure('Name','Tangent-Space Steps')

subplot(1,2,1)
stairs(continuationVector(:),[PC_steps(:);PC_steps(end)])
grid on
xlabel('Continuation Parameter','Interpreter','latex')
ylabel('Stepsize','Interpreter','latex')

subplot(1,2,2)
plot(continuationVector)
grid on
xlabel('PC Iterations','Interpreter','latex')
ylabel('Continuation Parameter','Interpreter','latex')
end