function plotControlParameters(xi,continuationVector)
figure('Name','Change of Control Parameters')
hold on
nXi = size(xi,2);

for iXi = 1:nXi
    plot(continuationVector,xi(:,iXi),'DisplayName',['$\xi_{',num2str(iXi),'}$'])
end
grid on
xlabel('Continuation Parameter','Interpreter','latex')
legend('Interpreter','latex')

end
