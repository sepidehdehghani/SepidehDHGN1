function[ETC,sol]=ObjectiveFunction(z,model)
% OBJECTIVEFUNCTION calculate the value of objective function according to
% decision variables and parameters

    global NFE;
    if isempty(NFE)
        NFE=0;
    end
    NFE=NFE+1;
%%%decision variables assignment
n=z(1);           %%%Sample size
h=z(2);           %%%the sampling intervals
% L=z(3);           %%%control chart limit
K=z(4);           %%% The number of samplings of monitoring period
% thetarma=z(5);    %%%  The autoregressive parameter of ARMA control chart
% phi=z(6);         %%% The  moving average parameter  of the ARMA control chart


%%% Parameters assignment
    nu=model.nu;
    cpm=model.cpm;
    ccm=model.ccm;
    cpm1=model.cpm1;
    ccm1=model.ccm1;
    crm1=model.crm1;
    Cy=model.Cy;
    ch=model.ch;
    cf=model.cf;
    cv=model.cv;
    Cd=model.Cd;
    CR=model.CR;
    P=model.P;
    Pr=model.Pr;
    d1=model.d1;
    D=model.D;
    s=model.s;
    lambda=model.lambda;
    lambda0=model.lambda0;
    landa1=model.landa1;
    landa2=model.landa2;
    lower=model.lower;
    upper=model.upper;
    theta0=model.theta0;
    theta1=model.theta1;
    gamma0=model.gamma0;
    gamma1=model.gamma1;
    mu0=model.mu0;
    sigma=model.sigma;
    mu1=model.mu1;
    USL=model.USL;
    LSL=model.LSL;
    T=(K+1).*h;

%% Time to shift variable
    f0=@(t)      (lambda0.*(nu.*t.^(nu-1)).*exp(-lambda0.*t.^nu));
    Fbar0=@(x)   exp(-lambda0.*x.^nu);

%%%middle parameters that depends on desision variables
%     alfa= 2*normcdf(-L);

%%%middle parameters that depend on the type of assignable cause
%     beta=zeros(1,s);
% for j=1:s
%     beta(j)=normcdf(L-delta(j)*sqrt(n))-normcdf(-L-delta(j)*sqrt(n));
% end

% %%%ARL
%     ARL0= 1/alfa;
%     ARL1= sum(lambda.*(1./(1-beta)))/lambda0;
% [ARL0,ARL1j]=SimulateArmaProcesswhile(z,model);
[ARL0,ARL1j]=SimulateArmaRandom(z,model);
ARL1= (sum(lambda.*ARL1j)/lambda0);

%% Time to failure variable in in-control state
    g0=@(x)       gamma0.*(theta0.*x.^(theta0-1)).*exp(-gamma0.*x.^theta0);
    Gbar0=@(t)    exp(-gamma0.*t.^theta0);
%% Scenario 1

%%% probability
    PS(1)=integral(@(x) g0(x).*Fbar0(x),0,T);

%%% Expected Time in control and out of control state
    ETinS(1)=integral(@(x) x.*g0(x).*Fbar0(x),0,T)./PS(1);
    EToutS(1)=0;
    ETPS(1)=ETinS(1)+EToutS(1);

%%% Expected non-conforming items
    Nin=1-(normcdf((USL-mu0)./sigma)-normcdf((LSL-mu0)./sigma));
    ENinS(1)=P.*Nin.*ETinS(1);
    ENoutS(1)=0;
    ENS(1)=ENinS(1)+ENoutS(1);

%%% Expected Rework cycle time
    ETrS(1)=(1-d1)*ENS(1)/Pr;

%%% Expected inventory time
    ECTS(1)=(P*ETPS(1)-ENS(1)+(1-d1)*ENS(1))/D;

%%% Expected cycle time for satisfy demand
    ETdS(1)=ECTS(1)-ETPS(1)-ETrS(1);

%%% Expected Inspection cost
    EICS(1)=(Cy.*ETPS(1))/(h.*ARL0);

%%% Expected Maintenance cost
    EMS(1)=ccm;

%%% Expected sampling cost
    ECSampleS(1)=(cf+n*cv)*ETPS(1)/h;

%%% Quality loss cost
    Cin=P*Nin*((1-d1)*CR+d1*Cd);
    EQCS(1)=Cin*ETPS(1);

%%% Expected inventory holding cost
    SQ11=((P-D)*ETPS(1)^2)/2;
    SQ12=(ETdS(1)*D+(P-D)*ETPS(1)-ENS(1))*ETrS(1)/2;
    SQ13=(D*ETdS(1)^2)/2;
    EHCS(1)=ch*(SQ11+SQ12+SQ13);

%%% Total expected cost
    ETCS(1)=EICS(1)+EMS(1)+ECSampleS(1)+EQCS(1)+EHCS(1);

%% Scenario 2
    PS(2)=Fbar0(T).*Gbar0(T);

%%% 0Expected Time in control and out of control state
    ETinS(2)=T;
    EToutS(2)=0;
    ETPS(2)=ETinS(2)+EToutS(2);

%%% Expected non-conforming items
    ENinS(2)=P.*Nin.*ETinS(2);
    ENoutS(2)=0;
    ENS(2)=ENinS(2)+ENoutS(2);

%%% Expected Rework cycle time
    ETrS(2)=(1-d1)*ENinS(2)/Pr;

%%% Expected inventory time
    ECTS(2)=(P*ETPS(2)-ENS(2)+(1-d1)*ENS(2))/D;

%%% Expected cycle time for satisfy demand
    ETdS(2)=ECTS(2)-ETPS(2)-ETrS(2);

%%% Expected Inspection cost
    EICS(2)=Cy.*K/ARL0;

%%% Expected Maintenance cost
    EMS(2)=cpm;

%%% Expected sampling cost
    ECSampleS(2)=(cf+n*cv)*K;

%%% Quality loss cost
    EQCS(2)=Cin*ETPS(2);

%%% Expected inventory holding cost
    SQ21=((P-D)*ETPS(2)^2)/2;
    SQ22=(ETdS(2)*D+(P-D)*ETPS(2)-ENS(2))*ETrS(2)/2;
    SQ23=(D*ETdS(2)^2)/2;
    EHCS(2)=ch*(SQ21+SQ22+SQ23);

%%% Total expected cost
    ETCS(2)=EICS(2)+EMS(2)+ECSampleS(2)+EQCS(2)+EHCS(2);

%% Scenario 3

%%%probability
    mult3=zeros(1,s);
    ETijS3=zeros(1,s);
    ETojS3=zeros(1,s);
    Nout=zeros(1,s);
    ENoS3=zeros(1,s);
    Cout=zeros(1,s);
    PS(3)=0;
 for j=1:s
    Gbar1j=@(t)   exp(-gamma1(j).*t.^theta1);
    g1j=@(y)      gamma1(j).*(theta1.*y.^(theta1-1)).*exp(-gamma1(j).*y.^theta1);
    mult3(j)=integral2(@(t,y) Gbar0(t).*f0(t).*((1-(1./ARL1j(j))).^ceil((y-t)/h)).*g1j(y)./Gbar1j(t),0,T,@(t) t,T,'AbsTol',1e-10,'RelTol',1e-1);
    PS(3)=PS(3)+(lambda(j).*mult3(j))./lambda0;

    ETijS3(j)=integral2(@(t,y) t.*Gbar0(t).*f0(t).*((1-(1./ARL1j(j))).^ceil((y-t)/h)).*g1j(y)./(Gbar1j(t).*mult3(j)),0,T,@(t) t,T,'AbsTol',1e-10,'RelTol',1e-1);
    ETojS3(j)=integral2(@(y,t) y.*Gbar0(t).*g1j(y).*((1-(1./ARL1j(j))).^ceil((y-t)/h)).*f0(t)./(Gbar1j(t).*mult3(j)),0,T,0,@(y) y,'AbsTol',1e-10,'RelTol',1e-1)-ETijS3(j);

    Nout(j)=1-(normcdf((USL-mu1(j))./sigma)-normcdf((LSL-mu1(j))./sigma));
    ENoS3(j)=Nout(j).*ETojS3(j);

    Cout(j)=P.*Nout(j)*((1-d1)*CR+d1*Cd);
 end
 
%%% Expected Time in control and out of control state
    ETinS(3)=sum(lambda.*ETijS3)/lambda0;
    EToutS(3)=sum(lambda.*ETojS3)/lambda0;
    ETPS(3)= ETinS(3)+EToutS(3);
   
%%% Expected non-conforming items
    ENinS(3)=P.*Nin.*ETinS(3);
    ENoutS(3)=P.*sum(lambda.*ENoS3)./lambda0;
    ENS(3)=ENinS(3)+ENoutS(3);
    
%%% Expected Rework cycle time
    ETrS(3)=(1-d1)*ENS(3)./Pr;

%%% Expected inventory time
    ECTS(3)=(P*ETPS(3)-ENS(3)+(1-d1)*ENS(3))/D;

%%% Expected cycle time for satisfy demand
    ETdS(3)=ECTS(3)-ETPS(3)-ETrS(3);

%%% Expected Inspection cost
    EICS(3)=(Cy.*ETPS(3))./(h.*ARL0);

%%% Expected Maintenance cost
    EMS(3)=(sum(lambda.*ccm1)/lambda0);

%%% Expected sampling cost
    ECSampleS(3)= (cf+n*cv)*ETPS(3)/h;

%%% Quality loss cost
    EQCS(3)=Cin*ETinS(3)+EToutS(3)*(sum(lambda.*Cout)/lambda0);

%%% Expected inventory holding cost
    SQ31=((P-D)*ETPS(3)^2)/2;
    SQ32=(ETdS(3)*D+(P-D)*ETPS(3)-ENS(3))*ETrS(3)/2;
    SQ33=(D*ETdS(3)^2)/2;
    EHCS(3)=ch*(SQ31+SQ32+SQ33);

%%% Total expected cost
    ETCS(3)=EICS(3)+EMS(3)+ECSampleS(3)+EQCS(3)+EHCS(3);
 

%% Scenario 4

%%%probability
    mult4=zeros(1,s);
    mult5=zeros(1,s);
    mult6=zeros(1,s);
    mult7=zeros(1,s);
    mult8=zeros(1,s);
    ETiS4=zeros(1,s);
    EToS4=zeros(1,s);
    Nout=zeros(1,s);
    ENoS4=zeros(1,s);
    PS(4)=0;
    for j=1:s 
    g1j=@(y)     gamma1(j).*(theta1.*y.^(theta1-1)).*exp(-gamma1(j).*y.^theta1);
    Gbar1j=@(t)  exp(-gamma1(j).*t.^theta1);
    mult4(j)=integral2(@(t,y) Gbar0(t).*f0(t).*g1j(y).*((1./ARL1j(j)).^ceil((y-t)./h))./Gbar1j(t),0,T,@(t) t,T,'AbsTol',1e-10,'RelTol',1e-1);
    mult5(j)=integral(@(t) Gbar0(t).*f0(t).*((1./ARL1j(j)).^ceil((T-t)./h)).*integral(@(y) g1j(y),T,inf)./Gbar1j(t),0,T);
    mult6(j)=(mult4(j)+mult5(j));
    PS(4)=PS(4)+(lambda(j).*mult6(j))./lambda0;

    mult7(j)=integral2(@(t,y) t.*Gbar0(t).*f0(t).*g1j(y).*((1./ARL1j(j)).^ceil((y-t)./h))./(Gbar1j(t).*mult6(j)),0,T,@(t) t,T,'AbsTol',1e-10,'RelTol',1e-1);
    mult8(j)=integral(@(t) t.*Gbar0(t).*f0(t).*((1./ARL1j(j)).^ceil((T-t)./h)).*integral(@(y) g1j(y),T,inf)./(Gbar1j(t).*mult6(j)),0,T);
    ETiS4(j)=mult7(j)+mult8(j);
    EToS4(j)=h.*ARL1j(j);
    
    Nout(j)=1-(normcdf((USL-mu1(j))./sigma)-normcdf((LSL-mu1(j))./sigma));
    ENoS4(j)=Nout(j).*EToS4(j);
   
    end
    
%%% Expected Time in control and out of control state
    ETinS(4)=sum(lambda.*ETiS4)./lambda0;
    EToutS(4)=sum(lambda.*EToS4)./lambda0;
    ETPS(4)=ETinS(4)+EToutS(4);
    
%%% Expected non-conforming items
    ENinS(4)=P.*Nin.*ETinS(4);
    ENoutS(4)=P.*sum(lambda.*ENoS4)./lambda0;
    ENS(4)=ENinS(4)+ENoutS(4);

%%% Expected Rework cycle time
    ETrS(4)=(1-d1)*ENS(4)/Pr;

%%% Expected inventory time
    ECTS(4)=(P*ETPS(4)-ENS(4)+(1-d1)*ENS(4))/D;

%%% Expected cycle time for fulfil demand
    ETdS(4)=ECTS(4)-ETPS(4)-ETrS(4);

%%% Expected Inspection cost
    EICS(4)=Cy.*(ETPS(4)/(h.*ARL0)+1);

%%% Expected Maintenance cost
    EMS(4)=(sum(lambda.*crm1))/lambda0;

%%% Expected sampling cost
    ECSampleS(4)=(cf+n*cv)*(ETinS(4)/h+sum(ARL1j.*lambda)/lambda0);

%%% Quality loss cost
    EQCS(4)=Cin*ETinS(4)+(sum(lambda.*Cout.*EToS4)/lambda0);

%%% Expected inventory holding cost
    SQ41=((P-D)*ETPS(4)^2)/2;
    SQ42=(ETdS(4)*D+((P-D)*ETPS(4)-ENS(4)))*ETrS(4)/2;
    SQ43=(D*ETdS(4)^2)/2;
    EHCS(4)=ch*(SQ41+SQ42+SQ43);

%%% Total expected cost
    ETCS(4)=EICS(4)+EMS(4)+ECSampleS(4)+EQCS(4)+EHCS(4);

%% Scenario5

%%% probability
    mult9=zeros(1,s);
    ETin5=zeros(1,s);
    EToS5=zeros(1,s);
    Nout=zeros(1,s);
    ENoS5=zeros(1,s);
    PS(5)=0;
    
for j=1:s
    g1j=@(y)  gamma1(j).*(theta1.*y.^(theta1-1)).*exp(-gamma1(j).*y.^theta1);
    Gbar1j=@(t)     exp(-gamma1(j).*t.^theta1);
    mult9(j)=integral(@(t) Gbar0(t).*f0(t).*((1-(1./ARL1j(j))).^ceil((T-t)./h)).*integral(@(y) g1j(y),T,inf)./Gbar1j(t),0,T);
    PS(5)=PS(5)+(lambda(j).*mult9(j))./lambda0;

    ETin5(j)=integral(@(t) t.*Gbar0(t).*f0(t).*((1-(1./ARL1j(j))).^ceil((T-t)./h)).*integral(@(y) g1j(y),T,inf)./(Gbar1j(t).*mult9(j)),0,T);
    EToS5(j)=T-ETin5(j);
    
    Nout(j)=1-(normcdf((USL-mu1(j))./sigma)-normcdf((LSL-mu1(j))./sigma));
    ENoS5(j)=Nout(j).*EToS5(j);

end

%%% Expected Time in control and out of control state
    ETinS(5)=sum(lambda.*ETin5)./lambda0;
    EToutS(5)=sum(lambda.*EToS5)./lambda0;
    ETPS(5)= ETinS(5)+EToutS(5);

%%% Expected non-conforming items
    ENinS(5)=P.*Nin.*ETinS(5);
    ENoutS(5)=P.*sum(lambda.*ENoS5)./lambda0;
    ENS(5)=ENinS(5)+ENoutS(5);

%%% Expected Rework cycle time
    ETrS(5)=(1-d1)*ENS(5)/Pr;

%%% Expected inventory time
    ECTS(5)=(P*ETPS(5)-ENS(5)+(1-d1)*ENS(5))/D;

%%% Expected cycle time for fulfil demand
    ETdS(5)=ECTS(5)-ETPS(5)-ETrS(5);

%%% Expected Inspection cost
    EICS(5)=Cy.*ETPS(5)/(h.*ARL0);

%%% Expected Maintenance cost
    EMS(5)=(sum(lambda.*cpm1)/lambda0);

%%% Expected sampling cost
    ECSampleS(5)= (cf+n*cv)*K;

%%% Quality loss cost
    EQCS(5)=Cin*ETinS(5)+EToutS(5)*(sum(lambda.*Cout)/lambda0);

%%% Expected inventory holding cost
    SQ51=((P-D)*ETPS(5)^2)/2;
    SQ52=(ETdS(5)*D+((P-D)*ETPS(5)-ENS(5)))*ETrS(5)/2;
    SQ53=(D*ETdS(5)^2)/2;
    EHCS(5)=ch*(SQ51+SQ52+SQ53);

%%% Total expected cost
    ETCS(5)=EICS(5)+EMS(5)+ECSampleS(5)+EQCS(5)+EHCS(5);
    ETIN=sum((ETinS.*PS)./ETPS);
    ETOUT=sum((EToutS.*PS)./ETPS);
    ENONC=sum(ENS.*PS);
    EIC=sum(EICS.*PS);
    EMC=sum(EMS.*PS);
    ESampleC=sum(ECSampleS.*PS);
    EQC=sum(EQCS.*PS);
    EHC=sum(EHCS.*PS);
% PS=PS1+PS2+PS3+PS4+PS5;
%% Mohasebe jarime 
violation1=max(0,1-(ARL0/lower));
violation2=max(0,(ARL1/upper)-1);

%%%%%%% Objective Function Calculation
ETC=ETCS(1)./ECTS(1)*PS(1)+ETCS(2)./ECTS(2)*PS(2)+ETCS(3)./ECTS(3)*PS(3)+ETCS(4)./ECTS(4)*PS(4)+ETCS(5)./ETCS(5)*PS(5)+landa1*violation1+landa2*violation2;

sol.PS=PS;
sol.ETPS=ETPS;
sol.ETinS=ETinS;
sol.EToutS=EToutS;
sol.ETC=ETC;
sol.EICS=EICS;
sol.EMS=EMS;
sol.ECSampleS=ECSampleS;
sol.EQCS=EQCS;
sol.EHCS=EHCS;
sol.ETIN=ETIN;
sol.ETOUT=ETOUT;
sol.ENONC=ENONC;
sol.T=T;
sol.ETrS=ETrS;
sol.ECTS= ECTS;
sol.z=z;
sol.ARL0=ARL0;
sol.ARL1=ARL1;
sol.EIC=EIC;
sol.EMC=EMC;
sol.ESampleC=ESampleC;
sol.EQC=EQC;
sol.EHC=EHC;
sol.violation1=violation1;
sol.violation2=violation2;
sol.IsFeasible=((violation1==0)&(violation2==0));
end
 