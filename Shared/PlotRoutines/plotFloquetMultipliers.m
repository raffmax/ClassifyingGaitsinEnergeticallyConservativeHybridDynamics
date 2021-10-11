function plotFloquetMultipliers(continuationVector,fmMatrix)
%plotFloquetMultipliers plots floquet multipliers in complex plane over
% a continuation parameter

if size(continuationVector,1) ~= size(fmMatrix,1)
    error('continuationVector must be the same number of columns as fmMatrix!')
end
nFM = size(fmMatrix,2); % number of floquet multiplier

%% plot 3D
figure('Name','Floquet Multipliers 3D');
hold on

[X,Y,Z] = cylinder(1,100);
Z = Z*(continuationVector(end)-continuationVector(1))+continuationVector(1);
surf(X,Y,Z,'EdgeColor','none','FaceColor',0.9*[1 1 1],'FaceAlpha',.5)
phi = linspace(0,2*pi,100); % for unit circle
plot3(sin(phi),cos(phi),continuationVector(1)*ones(1,100),'LineWidth',1.5,'color',0.3*[1 1 1])
plot3(sin(phi),cos(phi),continuationVector(end)*ones(1,100),'LineWidth',1.5,'color',0.3*[1 1 1])

for iFM = 1:nFM
    surface([real(fmMatrix(:,iFM))';real(fmMatrix(:,iFM))'],...
            [imag(fmMatrix(:,iFM))';imag(fmMatrix(:,iFM))'],...
            [continuationVector';continuationVector'],...
            [continuationVector';continuationVector'],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',4);  
end
colormap(flipud(winter))
caxis([continuationVector(1),continuationVector(end)])
colorbar

scatter3(real(fmMatrix(end,:)),...
         imag(fmMatrix(end,:)),...
         repmat(continuationVector(end),[1,nFM]),...
         150,'x','LineWidth',2,'MarkerEdgeColor',[0 0 0.6]); % reduced linewidth from 4 to 2

grid on
axis equal
xlabel('real','Interpreter','latex')
ylabel('imaginary','Interpreter','latex')
zlabel('Continuation Parameter','Interpreter','latex')

view(20,30)

end

