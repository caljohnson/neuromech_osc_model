function [Z2] = fsadj(Z,XLC)
% function - calculates time step for adjoint problem for HH model

  global gNa gK1 gK3 gL VNa VK VL cap dt
  
%  f(1)=-(gNa*XLC(2)^3*XLC(3)+gK3*XLC(4)^2+gK1*0.5^4+gL)/cap;
  f(1)=-(gNa*XLC(2)^3*XLC(3)+gK3*XLC(4)^2+gK1*XLC(5)^4+gL)/cap;
  f(2)=-(gNa*3*XLC(2)^2*XLC(3)*(XLC(1)-VNa))/cap;
  f(3)=-(gNa*XLC(2)^3*(XLC(1)-VNa))/cap;
  f(4)=-(gK3*2*XLC(4)*(XLC(1)-VK))/cap;
  f(5)=-(gK1*4*XLC(5)^3*(XLC(1)-VK))/cap;
  
  % ---- m ----
  c1=75.5; c2=1.0; c3=75.5; c4=13.5;
  % alpha = 40.0*(75.5.-1.0*X(1))/(exp(-(X(1)-75.5)/13.5)-1.0);
  alpha =40.0*(c1-c2*XLC(1))/(exp(-(XLC(1)-c3)/c4)-1.0);
  temp1 = exp(-(XLC(1)-c3)/c4)-1.0;
  ap=40.0*(-c2/temp1+(c1-c2*XLC(1))/temp1^2/c4*(temp1+1.0));
  
  c1=1.2262; c2=42.248;
  % beta = 1.2262/exp(X(1)/42.248);
  beta = c1/exp(XLC(1)/c2);
  bp=-c1/exp(XLC(1)/c2)/c2;

  g(2)=-(alpha+beta);
  h(2)= ap*(1.0-XLC(2))-bp*XLC(2);

  % ---- h ----
  c1=0.0035; c2=24.186;
  % alpha = 0.0035/exp(X(1)/24.186);
  alpha = c1/exp(XLC(1)/c2);
  ap=-c1/exp(XLC(1)/c2)/c2;

  c1=-51.25; c2=1.0; c3=-51.25; c4=5.2;
  % beta =  0.017*(-51.25-X(1))/(exp(-(51.25+X(1))/5.2)-1.0);
  temp1 = exp(-(XLC(1)-c3)/c4)-1.0;
  beta =0.017*(c1-c2*XLC(1))/temp1;
  bp=0.017*(-c2/temp1+(c1-c2*XLC(1))/temp1^2/c4*(temp1+1.0));
 
  g(3)=-(alpha+beta);
  h(3)= ap*(1.0-XLC(3))-bp*XLC(3);

  % ---- n3 ----
  c1=95.0; c2=1.0; c3=95.0; c4=11.8;
  %  alpha = (95.-X(1))/(exp(-(X(1)-95.)/11.8)-1.);
  temp1 = exp(-(XLC(1)-c3)/c4)-1.0;
  alpha =(c1-c2*XLC(1))/temp1;
  ap=-c2/temp1+(c1-c2*XLC(1))/temp1^2/c4*(temp1+1.0);
  
  c1=0.025; c2=22.222;
  % beta = 0.025/exp(X(1)/22.222);
  beta = c1/exp(XLC(1)/c2);
  bp=-c1/exp(XLC(1)/c2)/c2;

  g(4)=-(alpha+beta);
  h(4)= ap*(1.0-XLC(4))-bp*XLC(4);

  % ---- n1 ----
  c1=-44.0; c2=1.0; c3=-44.0; c4=2.3;  
  % alpha = 0.014*(-44.0-1.0*X(1))/(exp(-(X(1)-(-44.0))/2.3)-1.0);
  temp1 = exp(-(XLC(1)-c3)/c4)-1.0;
  alpha =0.014*(c1-c2*XLC(1))/temp1;
  ap=0.014*(-c2/temp1+(c1-c2*XLC(1))/temp1^2/c4*(temp1+1.0));
  
  c1=0.0043; c2=34.0; c3=44.0;
  % beta = 0.0043/exp((X(1)+44.0)/34.0); 
  beta = c1/exp((XLC(1)+c3)/c2);
  bp=-c1/exp((XLC(1)+c3)/c2)/c2;

  g(5)=-(alpha+beta);
  h(5)= ap*(1.0-XLC(5))-bp*XLC(5);


  dZ(1)=-(f(1)*Z(1)+h(2)*Z(2)+h(3)*Z(3)+h(4)*Z(4)+h(5)*Z(5));
  for i=2:5
%  dZ(1)=-(f(1)*Z(1)+h(2)*Z(2)+h(3)*Z(3)+h(4)*Z(4));
%  for i=2:4
      dZ(i)=-(f(i)*Z(1)+g(i)*Z(i));
  end;
  Z2(1:5)=Z(1:5)+dZ(1:5)*(-dt);
%    Z2(1:4)=Z(1:4)+dZ(1:4)*(-dt);

 
  
      