function[ARL0,ARL1j]=SimulateArmaRandom(z,model)
% Simulate ARMA Process
L=z(3);           %%%control chart limit
theta=z(5);    %%%  The autoregressive parameter of ARMA control chart
phi=z(6);         %%% The  moving average parameter  of the ARMA control chart

%%% Parameters assignment
    u=model.u;
    v=model.v;
    sigmaa=model.sigmaa;
    sigma=model.sigma;
    mu1=model.mu1;
    s=model.s;
    q=500;
mu0=model.mu0;
sigmaz=sqrt(((2*(theta-phi)*(1+theta)/(1+phi))+1)*sigma);
UCL=mu0+L*sigmaz;
LCL=mu0-L*sigmaz;
theta0=1+theta-phi;
r=10;
%%%%%%%%%%%%%%%%%
%%ARL0
Z=zeros(r,1);
X=zeros(r,1);

for l=1:q
    l
RL0(l)=0;
e=-inf;
while e<=0
      X=normrnd(mu0,sigma);
%     X(1)=normrnd(mu0,sigma);
%     w=normrnd(0,sigmaa,[r 1]);
%     Z(1)=mu0;
%     for i=2:r
%         X(i)=w(i)+v.*w(i-1)+u.*X(i-1);
%         Z(i)=theta0.*X(i)+theta.*X(i-1)+phi.*Z(i-1);
%     end 
%     a=max(0,Z(r)/UCL-1);
%     b=max(0,1-Z(r)/LCL);
    a=max(0,X/UCL-1);
    b=max(0,1-X/LCL);
    d=a+b;
    e=log(1+d);
    % d=1-xor(a,b);
    RL0(l)=RL0(l)+1;
end 
 end
%%ARL1
Z1=zeros(r,s);
Xp=zeros(r,1);
for j=1:s
    for l=1:q
    RL1(l,j)=0;
    o=-inf;
    while o<=0
    Xp(1,j)=normrnd(mu1(j),sigma);
    wp=normrnd(0,sigmaa,[r 1]);
    Z1(1,j)=mu1(j);
    for i=2:r
    Xp(i,j)=wp(i)+v.*wp(i-1)+u.*Xp(i-1,j);
    Z1(i,j)=theta0.*Xp(i,j)-theta.*Xp(i-1,j)+phi.*Z1(i-1,j);    
    end 
    f=max(0,Z1(r,j)/UCL-1);
    g=max(0,1-Z1(r,j)/LCL);
    m=f+g;
    o=log(1+m);
    % d=1-xor(a,b);
    RL1(l,j)=RL1(l,j)+1;
    end   
    end  
end
ARL0=mean(RL0);
ARL1j=mean(RL1);
end


