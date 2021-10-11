function createAnimationContinuationHopper(dataIN,type,fixedParam,full3D,boundaryRep)
%createAnimationContinuationPogoStick animates continuation and creates
%video
%   dataIN:          output of LC.getContinuationTrajectories
%   type:            'epsilon', 'E0'
%   fixedParam:      symbol, value, unit of fixed parameter in continuation
%   full3D:          logical value to create animation in 3D, otherwise 2D
%   boundaryRep:     repitions of limit cycles at start and end of continuation


if full3D % change the 3D lighting effects
    edges         = false; % show edges for torso and legs
else
    edges         = true; % show edges for torso and legs
end

% Colours [R,G,B] for torso;spring;foot
colours = [230, 230, 230; 255, 0, 0; 0, 0, 0]/255; 

% time-steps
dt = 7e-2;

% check if we simulate with or without pitch
if size(dataIN(1,1).x,2) == 8
   % add states
   % stateIn  = [x,y,alpha,l]
   % stateOut = [x,y,phi,alpha,l]
   data = dataIN;
   for iEps = 1:size(data,1)
       for iDomain = 1:size(data,2)
           n    = length(data(iEps,iDomain).t);
           col0 = zeros(n,1);
           data(iEps,iDomain).x  = [dataIN(iEps,iDomain).x(:,1:2),col0,dataIN(iEps,iDomain).x(:,3:4),...
                                     dataIN(iEps,iDomain).x(:,5:6),col0,dataIN(iEps,iDomain).x(:,7:8)];
       end
   end
else
    data = dataIN;
end

epsilon = 0;
hopper_graphic = HoppergraphicX(full3D,edges,colours,epsilon);
fig = gcf; % get figure
fig.Position =  [1 1 1000 1000];

% check screen size
screensize = get( groot, 'Screensize' );
if min(screensize(3:4)) < 1000
    warning('Your screen size is too small to cover the required 1000x1000 pixels!') 
end

nDomains   = size(data,2); % number of domains
nCon       = size(data,1); % number of continuation steps
F_complete = []; % figure data
contIDX    = [1*ones(1,nDomains*boundaryRep),...
              1:nCon,...
              nCon*ones(1,nDomains*boundaryRep)];

iDomain   = 0;
loopCount = 0;
for iCon = contIDX
    loopCount   = loopCount+1;
    iDomain     = iDomain+1;
    if iDomain > nDomains
        iDomain = 1; % reset
    end
    
    % z = motion(:, iDomain);
    tData = data(iCon,iDomain).t;
    t     = tData(1):dt:tData(end-1);
    x     = interp1(tData,data(iCon,iDomain).x,t);
    
    switch type
        case  'epsilon'
            if loopCount == length(contIDX)
                epsilonArray = linspace(data(iCon,1).epsilon,data(iCon,1).epsilon,length(t));
            else
                epsilonArray = linspace(data(iCon,1).epsilon,data(contIDX(loopCount+1),1).epsilon,length(t));
            end
            E0Array = zeros(1,length(t));
        case {'E0','E_avg','xDot_avg'}
            if loopCount == length(contIDX)
                E0Array      = linspace(data(iCon,1).(type),data(iCon,1).(type),length(t));
            else
                E0Array      = linspace(data(iCon,1).(type),data(contIDX(loopCount+1),1).(type),length(t));
            end
            epsilonArray = fixedParam.value*ones(1,length(t));
    end
    
    
    F(length(t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    
    for i=1:length(t)
        state_  = x(i,:);
        epsilon = epsilonArray(i);
        E0      = E0Array(i);

        update(hopper_graphic, state_', epsilon);

        switch type
            case 'epsilon'
                txtEps   = text(state_(1)-1.2,0,2.2,['$\varepsilon= ',num2str(round(epsilon,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                txtParam = text(state_(1)+.15,0,2.2,['$',fixedParam.symbol,' = ',num2str(round(fixedParam.value,2)),'~',fixedParam.unit,'$'],'FontSize',36,'Color','black','Interpreter','latex');
                stretchRectanglePlot = false;
            case {'E0','E_avg','xDot_avg'}
                switch type
                    case 'E0'
                        txtParam1 = text(state_(1)-1.2,0,2.2,['$E_0 = ',num2str(round(E0,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                    case 'xDot_avg'
                        txtParam1 = text(state_(1)-1.2,0,2.2,['$\dot{x}_\mathrm{avg} = ',num2str(round(E0,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                    otherwise
                        txtParam1 = text(state_(1)-1.2,0,2.2,['$E_\mathrm{avg} = ',num2str(round(E0,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                end
                % 
                if strcmp(type,'xDot_avg')
                    txtParam2 = text(state_(1)-0.55,0,2.2,'$\sqrt{g l_o}$','FontSize',36,'Color','blue','Interpreter','latex');
                else
                    txtParam2 = text(state_(1)-0.55,0,2.2,'$m g l_o$','FontSize',36,'Color','blue','Interpreter','latex');
                end
                txtEps    = text(state_(1)+.15,0,2.2,['$\varepsilon= ',num2str(round(epsilon,2)),'$'],'FontSize',36,'Color','black','Interpreter','latex');
                stretchRectanglePlot = true;
        end
        txtT   = text(state_(1)+.15,0,2.0,['$T~= ',sprintf('%.2f',data(iCon,end).t(end)),'~\sqrt{l_o/g}$'],'FontSize',36,'Color','black','Interpreter','latex');
        % plot rectangle
        x1 = state_(1)-1.25;
        if ~stretchRectanglePlot
            x2 = state_(1)-0.65;
        else
            x2 = state_(1)-0.25;
        end
        x3 = x2;
        x4 = x1;
        y  = 0;
        z1 = 2.3;
        z2 = 2.3;
        z3 = 2.1;
        z4 = 2.1;
        hold on
        rectangle = plot3( [x1 x2 x3 x4 x1], [y y y y y], [z1 z2 z3 z4 z1],'b');

        % append
        F(frame_counter) = getframe(fig);
        if stretchRectanglePlot
            delete(txtEps)
            delete(txtT)
            delete(txtParam1)
            delete(txtParam2)
            delete(rectangle)
        else
            delete(txtEps)
            delete(txtT)
            delete(txtParam)
            delete(rectangle)
        end
        frame_counter = frame_counter + 1;

    end
    F_complete = [F_complete, F];
    clearvars F % s.t. F can be resized

end

video_gif_ID = '\pogoStick';
v = VideoWriter(strcat(pwd, video_gif_ID, '.mp4'),'MPEG-4');
open(v)
writeVideo(v,F_complete)
close(v)
clearvars v 
close gcf
end