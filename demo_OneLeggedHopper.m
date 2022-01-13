% DEMO_ONELEGGEDHOPPER A demo on how to use ModelHomotopy.
%   This script demonstrates how to generate an object of OneLeggedHopper 
%   which is of type MechanicalModel and how to set up a root finding 
%   problem for an isolated solution (limit cycle).
%   It uses pseudo-arclength continuation in the exploration of connected
%   strata.
%
%   See also Solution, OneLeggedHopper, RootFindingProblem, LimitCycle,
%            PseudoArclengthPC
clear, clc, close all;
install_project();
%% load model
% - sub-class of MechanicalModel -
hopper = OneLeggedHopper('load', false);

%% initial guess
% q = [x,y,phi,alpha,l]
yApex = 1.001;

% get model parameters to compute domain times
m   = hopper.parameters.m;
g   = hopper.parameters.g;
k   = hopper.parameters.param_k_l(0);
l0  = hopper.parameters.l_0;

% flight phase
tF = sqrt(2*(yApex-l0)/g); 

% stance phase: y = a*sin( w*t + phi ) + l0 - g/(w^2)
w   = sqrt(k/m);
a   = g/w^2*sqrt(tF^2*w^2+1);
phi = acos(-g*tF/(a*w));
tS  = (acos(g*tF/(a*w))-phi+2*pi)/w;

x_init = [0; l0; 0;0; l0; 0; -g*tF; 0;0; -g*tF]; 
% domain times: flight, stance, flight
t_init  = [tS, 2*tF]';
% control parameter (numerical damping)
xi_init = 0;
% model homotopy parameter
epsilon = 0;
% store initial results in Solution
sol = Solution(x_init,t_init,xi_init,epsilon);
sol.param.E0 = yApex;

% configure sequence:
% for this model, we have a constraint during flight and a constraint
% during stance, thus the phase variable is dim(z) = 2
% z = [stance constraint active;....
%      flight constraint active]
sequence.FlowMapOrder   = [1, 0;... % stance constraint
                           0, 1];   % flight constraint
% the jump map follows after the flowMap and needs the active
% constraint configuration of the next flowMap.
% Thus, JumpMapOrder is usually just a shifted FlowMapOrder.
sequence.JumpMapOrder   = [0, 1;...
                           1, 0];
% define sequence of events, that should be triggered at domain times t_i
% the specification (here 1,2) is user-defined and corresponds to the event
% implementation in their MechanicalModel
sequence.events         = [2, 1]; % lift-off, touch-down

% fixing the state to its inital value is done here.
fixedStates             = [1,3,8]; % fixed at initial state
% define the states which need to be periodic
periodicStates          = [2,4,5,6,7,9,10];

% set anallytic first order gradients
opts.Grad1 = true; 

% set up RootSearch to define the search for a limit cycle
rfp = RootFindingProblem(hopper,sequence,'fixedStates',fixedStates,...
                                        'periodicStates',periodicStates,...
                                        'simulation',Simulation('odeSolver',@ode113),...
                                        'options',opts);
% create limit cycle
solAn  = rfp.solveNewtonRFP(sol);
LC     = LimitCycle(['E0=',num2str(yApex)],solAn,rfp);

%% run Continuation in E0= E bar
tic
% jump over simple bifurcations to explore M_0, M_1 and M_4
arcPC  = PseudoArclengthPC(rfp,LC);
target = 2.4; % E_bar = E_0 = 2.4
[exitflag,PC_steps] = arcPC.runContinuation('E0',target);

% approximate 2 simple bifurcations
[LCBefore,LCAfter] = LC.findNextBifurcation();
[LCBif1, flag1]    = arcPC.approximateRoot(LCBefore,LCAfter, 'detJacAug'); 
[LCBefore,LCAfter] = LCBif1.findNextBifurcation();
[LCBif2, flag2]    = arcPC.approximateRoot(LCBefore,LCAfter, 'detJacAug'); 

% get new gaits on new generators
[LC_M2,LC_M3] = arcPC.branchingOff(LCBif1);
[LC_M5,LC_M6] = arcPC.branchingOff(LCBif2);

% initialize new generators
M2 = copy(LC_M2);
M2.disconnectLC;
M3 = copy(LC_M3);
M3.disconnectLC;
M5 = copy(LC_M5);
M5.disconnectLC;
M6 = copy(LC_M6);
M6.disconnectLC;

% explore generators
% save trajectory of gaits at E_0 = E_bar = 1.8
GaitTrajectories_18 = cell(7,1);
M = {M2,M3,M5,M6};
l = [2 , 3, 5, 6];
for i = 1:length(M)
    arcPC_Mi  = PseudoArclengthPC(rfp,M{i});
    target = 1.8; % E_bar = E_0 = 1.8
    arcPC_Mi.runContinuation('E0',target);
    M{i} = M{i}.getLCend;
    % get trajectory
    GaitTrajectories_18{l(i)} = M{i}.getTrajectory(3e-2);
    arcPC_Mi  = PseudoArclengthPC(rfp,M{i});
    target = 2.4; % E_bar = E_0 = 2.4
    arcPC_Mi.runContinuation('E0',target); 
end

% find E_0 = Ebar = 1.8 on M4
lc4 = LimitCycle('',arcPC.Bifurcations(2).Next(1).Sol,rfp);
arcPC_4  = PseudoArclengthPC(rfp,lc4);
target = 1.8; % E_bar = E_0 = 1.8
arcPC_4.runContinuation('E0',target);
GaitTrajectories_18{4} = lc4.getTrajectory(3e-2);

toc

%% extract data from strata
M0 = LC.getContinuationData;
M1 = arcPC.Bifurcations(1).Next(1).getContinuationData;
M2 = M2.getLCstart.getContinuationData;
M3 = M3.getLCstart.getContinuationData;
M4 = arcPC.Bifurcations(2).Next(1).getContinuationData;
M5 = M5.getLCstart.getContinuationData;
M6 = M6.getLCstart.getContinuationData;

BP1 = arcPC.Bifurcations(1);
BP2 = arcPC.Bifurcations(2);

clc; % clear command window

%% Plots
figure
grid on
hold on
plot([M0.E0;BP1.Sol.param.E0],[sum(M0.tDomain,2);sum(BP1.Sol.tDomain)],'b')
plot([BP1.Sol.param.E0;M1.E0;BP2.Sol.param.E0],[sum(BP1.Sol.tDomain);sum(M1.tDomain,2);sum(BP2.Sol.tDomain)],'g')
plot([BP1.Sol.param.E0;M2.E0],[sum(BP1.Sol.tDomain);sum(M2.tDomain,2)],'r')
plot([BP1.Sol.param.E0;M3.E0],[sum(BP1.Sol.tDomain);sum(M3.tDomain,2)],'r--')
plot([BP2.Sol.param.E0;M4.E0],[sum(BP2.Sol.tDomain);sum(M4.tDomain,2)],'c')
plot([BP2.Sol.param.E0;M5.E0],[sum(BP2.Sol.tDomain);sum(M5.tDomain,2)],'k')
plot([BP2.Sol.param.E0;M6.E0],[sum(BP2.Sol.tDomain);sum(M6.tDomain,2)],'k--')
scatter(BP1.Sol.param.E0,sum(BP1.Sol.tDomain),'k','filled')
scatter(BP2.Sol.param.E0,sum(BP2.Sol.tDomain),'k','filled')
xlabel('$\bar{E}$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
legend('$\mathcal{M}_0$','$\mathcal{M}_1$','$\mathcal{M}_2$','$\mathcal{M}_3$',...
       '$\mathcal{M}_4$','$\mathcal{M}_5$','$\mathcal{M}_6$','BC1','BC2','Interpreter','latex')
                           
figure
grid on
hold on
plot([M0.E0;BP1.Sol.param.E0],[M0.x0(:,6);BP1.Sol.x0(6)],'b')
plot([BP1.Sol.param.E0;M1.E0;BP2.Sol.param.E0],[BP1.Sol.x0(6);M1.x0(:,6);BP2.Sol.x0(6)],'g')
plot([BP1.Sol.param.E0;M2.E0],[BP1.Sol.x0(6);M2.x0(:,6)],'r')
plot([BP1.Sol.param.E0;M3.E0],[BP1.Sol.x0(6);M3.x0(:,6)],'r--')
plot([BP2.Sol.param.E0;M4.E0],[BP2.Sol.x0(6);M4.x0(:,6)],'c')
plot([BP2.Sol.param.E0;M5.E0],[BP2.Sol.x0(6);M5.x0(:,6)],'k')
plot([BP2.Sol.param.E0;M6.E0],[BP2.Sol.x0(6);M6.x0(:,6)],'k--')
xlabel('$\bar{E}$','Interpreter','latex')
ylabel('$\dot{x}_0$','Interpreter','latex')
scatter(BP1.Sol.param.E0,BP1.Sol.x0(6),'k','filled')
scatter(BP2.Sol.param.E0,BP2.Sol.x0(6),'k','filled')
legend('$\mathcal{M}_0$','$\mathcal{M}_1$','$\mathcal{M}_2$','$\mathcal{M}_3$',...
       '$\mathcal{M}_4$','$\mathcal{M}_5$','$\mathcal{M}_6$','BC1','BC2','Interpreter','latex')
%% create animations for gaits on different strata at E_bar = 1.8
% uses animation from https://bitbucket.org/ganzhenyu/bipedalslipmodel_swingleg
framesM2 = getAnimationHopper(GaitTrajectories_18{2});
framesM3 = getAnimationHopper(GaitTrajectories_18{3});
framesM4 = getAnimationHopper(GaitTrajectories_18{4});
framesM5 = getAnimationHopper(GaitTrajectories_18{5});
framesM6 = getAnimationHopper(GaitTrajectories_18{6});

% export as gif
movie2gif(framesM2, 'M2_OneLeggedHopper.gif','DelayTime',.05,'LoopCount',Inf)
movie2gif(framesM3, 'M3_OneLeggedHopper.gif','DelayTime',.05,'LoopCount',Inf)
movie2gif(framesM4, 'M4_OneLeggedHopper.gif','DelayTime',.05,'LoopCount',Inf)
movie2gif(framesM5, 'M5_OneLeggedHopper.gif','DelayTime',.05,'LoopCount',Inf)
movie2gif(framesM6, 'M6_OneLeggedHopper.gif','DelayTime',.05,'LoopCount',Inf)
