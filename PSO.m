clc;
clear;
close all;

%% Problem Definition

global NFE;
model=CreateModelFirst();

CostFunction=@(z) ObjectiveFunction(z,model); % Cost Function


nVar=6;                          % Number of Decision Variables
VarSize=[1 nVar];                % Size of Decision Variables Matrix

VarMin=[0 0.4 2.5 0 .01 .001];              % Lower Bound of Variables
VarMax=[1 2 4.5 1 1 1];                   % Upper Bound of Variables


%% PSO Parameters

MaxIt=1;         % Maximum Number of Iterations

nPop=20;           % Population Size (Swarm Size)

% w=1; 
% % Inertia Weight
% wdamp=0.5;     % Inertia Weight Damping Ratio
% c1=2;           % Personal Learning Coefficient
% c2=2;           % Global Learning Coefficient

% Constriction Coefficients
phi1=2.05;
phi2=2.05;
phi=phi1+phi2;
chi=2/(phi-2+sqrt(phi^2-4*phi));
w=chi;          % Inertia Weight
wdamp=1;        % Inertia Weight Damping Ratio
c1=chi*phi1;    % Personal Learning Coefficient
c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

%% Initialization

empty_particle.Position=[];
empty_particle.PositionInt=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.sol=[];
empty_particle.Best.Position=[];
empty_particle.Best.PositionInt=[];
empty_particle.Best.Cost=[];
empty_particle.Best.sol=[];


particle=repmat(empty_particle,nPop,1);

GlobalBest.Cost=inf;

for i=1:nPop
    
    % Initialize Position
    particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    particle(i).PositionInt(1)=min(floor(20*particle(i).Position(1)+1),20);
    particle(i).PositionInt(2)=particle(i).Position(2);
    particle(i).PositionInt(3)=particle(i).Position(3);
    particle(i).PositionInt(4)=min(floor(39+31*particle(i).Position(4)+1),70);
    particle(i).PositionInt(5)=particle(i).Position(5);
    particle(i).PositionInt(6)=particle(i).Position(6);
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    [particle(i).Cost, particle(i).sol]=CostFunction(particle(i).PositionInt);
    
    % Update Personal Best
    particle(i).Best.Position=particle(i).Position;
    particle(i).Best.PositionInt=particle(i).PositionInt;
    particle(i).Best.Cost=particle(i).Cost;
    particle(i).Best.sol=particle(i).sol;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end

BestCost=zeros(MaxIt,1);

nfe=zeros(MaxIt,1);


%% PSO Main Loop
ETTC=zeros(1,MaxIt);
for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        
       particle(i).PositionInt(1)=min(floor(20*particle(i).Position(1)+1),20);
       particle(i).PositionInt(2)=particle(i).Position(2);
       particle(i).PositionInt(3)=particle(i).Position(3);
       particle(i).PositionInt(4)=min(floor(39+31*particle(i).Position(4)+1),70);
       particle(i).PositionInt(5)=particle(i).Position(5);
       particle(i).PositionInt(6)=particle(i).Position(6);
        % Evaluation
        [particle(i).Cost, particle(i).sol]=CostFunction(particle(i).PositionInt);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.PositionInt=particle(i).PositionInt;
            particle(i).Best.Cost=particle(i).Cost;
            particle(i).best.sol= particle(i).sol;
            
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    ETTC(it)=GlobalBest.Cost;
    nfe(it)=NFE;
    disp(['iter=  ' num2str(it) ': NFE = ' num2str(nfe(it)) ',ETC=  ' num2str(GlobalBest.Cost)])
    disp(GlobalBest.PositionInt)
    
    w=w*wdamp;
    
end

%% Results
itter=1:MaxIt;
 plot(itter,ETTC,'r');
%  figure;
% plot(BestCost);
% % semilogy(nfe,GlobalBest.Cost,'LineWidth',2);
% xlabel('Iteration');
% ylabel('GlobalBest.Cost');

% n=20;
% L=1:4;
% s=7;
% delta=0.25;
% % deltaV= 0.09;
% % delta=deltaB:deltaV:deltaB+(s-1)*deltaV;

% lambda=0.0016;
% % lambda=[0.0008 0.0004 0.0002 0.0001 0.00005 0.000025 0.0000125]; %scale parameter for assignable cause
% % lambda0=sum(lambda);
%     alfa= 2*normcdf(-L);
% 
% %%%middle parameters that depend on the type of assignable cause
% %     beta=zeros(1,s);
% 
%     beta=normcdf(L-delta*sqrt(n))-normcdf(-L-delta*sqrt(n));
% 
% 
% %%ARL
%     ARL0= 1./alfa;
%     ARL1=1./(1-beta);
% lower=200;
% upper=5;
% % L=1:4;
% y=ARL0/lower;
% % y>1;
% yp=ARL1/upper;
% % yp<1;
% plot(L,y,'r');
% hold on;
% plot(L,yp,'b');
% 
% 
% 
