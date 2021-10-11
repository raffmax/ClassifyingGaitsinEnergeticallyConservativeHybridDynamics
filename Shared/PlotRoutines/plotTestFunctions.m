function plotTestFunctions(conData,continuationVector)
figure('Name','Test-Functions for Bifurcations')
counter = 1;
if isfield(conData,'detHessian')
    optimize = true;
    nRow     = 3;
    nCol     = 3;
else
    optimize = false;
    nRow     = 2;
    nCol     = 3;
end
    
subplot(nRow,nCol,counter)
plot(continuationVector,conData.tDomain)
grid on
xlabel('Continuation Parameter','Interpreter','latex')
title('Domain Times')
counter = counter + 1;

if ~isempty(conData.eventVelocities)
    subplot(nRow,nCol,counter)
    plot(continuationVector,conData.eventVelocities)
    grid on
    xlabel('Continuation Parameter','Interpreter','latex')
    title('Event Grazing')
    counter = counter + 1;
end

subplot(nRow,nCol,counter)
plot(continuationVector,conData.detSimpleJacobi)
grid on
xlabel('Continuation Parameter','Interpreter','latex')
title('Simple Jacobian')
counter = counter + 1;

subplot(nRow,nCol,counter)
plot(continuationVector,conData.detAugJacobi)
grid on
xlabel('Continuation Parameter','Interpreter','latex')
title('Augmented Jacobian')
counter = counter + 1;

% if isfield(conData,'FloquetMultipliers')
%     subplot(nRow,nCol,counter)
%     plot(continuationVector,abs(conData.FloquetMultipliers(:,2:end))-1)
%     grid on
%     xlabel('Continuation Parameter','Interpreter','latex')
%     title('Absolute Floquet Multipliers')
%     counter = counter + 1;
% end

if optimize
    subplot(nRow,nCol,counter)
    plot(continuationVector,conData.detConJacobi)
    grid on
    xlabel('Continuation Parameter','Interpreter','latex')
    title('Regularity of Equality Constraints')
    counter = counter + 1;
    
    subplot(nRow,nCol,counter)
    plot(continuationVector,min(real(conData.suffCondEig),[],2))
    grid on
    xlabel('Continuation Parameter','Interpreter','latex')
    title('Sufficient Condition')
    counter = counter + 1;
end

end