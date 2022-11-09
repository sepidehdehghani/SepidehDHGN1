function model=CreateModelFirst()
% the values of parameters
s=7;                      % number of assignable cause
nu=1.7;                   % shape parameter of assignable cause

d1=0.1;                   % rework                                   
cv=0.1;                   % variable cost of sampling
ch=10;                    % inventory holding cost per unit 
ccm=3900;                 % CM cost for scenario 1
ccmB=2000;               % CM cost for scenario 3,5
ccV= 500;
ccm1= ccmB:ccV:ccmB+(s-1)*ccV;

cf=0.5;                   % fix cost of samoling
Cy= 400;                  % inspection cost per unit
cpm= 1300;                % PM cost for Scenario 2
mu0=5;                  % The process mean in-control state
theta0=1.4;               % shape parameter of failure in control state 
theta1=1.4;               % shape parameter of failure out of control state
gamma0=0.0001;            % scale parameter of failure in control state 
Cd=20;                    % scrape  production cost per unit
CR=5;                     % rework Cost per unit 
P=100;                    % production rate 
Pr=150;                   % rework rate 
d=10000;                  % The annual demand rate
D=50;                     % Demand rate in units per unit of time
sigma=1.08;                  % The standard deviation of quality characteristic
USL=3;                    % upper bound of The process mean
LSL=1;                    % lower bound of The process mean

deltaB=0.56;              % The magnitude of the shift in the mean of quality characteristic 
deltaV= 0.03;
delta=deltaB:deltaV:deltaB+(s-1)*deltaV;

lambdaB=0.1;             % Scale parameter of the Weibull distribution assignable cause
pow=1:s;
lambda= lambdaB./(2.^(pow-1));
lambda0=sum(lambda);

gammaB=0.14;              % Scale parameter of the Weibull distribution failure in out of control
pow=1:s;
gamma1=gammaB./(2.^(pow-1));

mu1=mu0+sigma.*delta;     % The process mean out of control state

cpmB= 1500;               % PM cost for Scenario 5
cpV= 400;
cpm1= cpmB:cpV:cpmB+(s-1)*cpV;

crmB= 1200;               % PRDM cost 
crV= 300;
crm1=crmB:crV:crmB+(s-1)*crV;


landa1=100;
landa2=100;

lower=300;
upper=10;

%%% Model assignment
model.nu=nu;
model.cpm=cpm;
model.ccm=ccm;
model.cpm1=cpm1;
model.ccm1=ccm1;
model.crm1=crm1;
model.Cy=Cy;
model.ch=ch;
model.cf=cf;
model.cv=cv;
model.Cd=Cd;
model.CR=CR;
model.P=P;
model.Pr=Pr;
model.d=d;
model.d1=d1;
model.D=D;
model.s=s;
model.lambda=lambda;
model.delta= delta;
model.lambda0= lambda0;
model.landa1=landa1;
model.landa2=landa2;
model.lower=lower;
model.upper=upper;
model.gamma0=gamma0;
model.gamma1=gamma1;
model.theta0=theta0;
model.theta1=theta1;
model.mu0=mu0;
model.sigma=sigma;
model.mu1=mu1;
model.USL=USL;
model.LSL=LSL;
end