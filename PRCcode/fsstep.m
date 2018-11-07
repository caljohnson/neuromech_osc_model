function [Y] = fsstep(X,cbias)
% function - calculates time step for FS model

  global gNa gK1 gK3 gL VNa VK VL cap dt

  INa = gNa*X(2)^3*X(3)*(X(1)-VNa);
  IK = (gK3*X(4)^2+gK1*X(5)^4)*(X(1)-VK);
  IL = gL*(X(1)-VL);
  Y(1) = X(1) + (cbias - IL - INa - IK)*dt/cap;  

  % -- m --
  alpha = 40.0*(75.5-X(1))/(exp(-(X(1)-75.5)/13.5)-1.0);
  beta = 1.2262/exp(X(1)/42.248);
  Y(2) = X(2) + (alpha*(1.0-X(2)) - beta*X(2))*dt;
  
  % -- h --
  alpha = 0.0035/exp(X(1)/24.186);
  beta =  -0.017*(51.25+X(1))/(exp(-(51.25+X(1))/5.2)-1.0);
  Y(3) = X(3) + (alpha*(1.0-X(3)) - beta*X(3))*dt;
  
  % -- n3 --
  alpha = (95.-X(1))/(exp(-(-95.+X(1))/11.8)-1.);
  beta = 0.025/exp(X(1)/22.222);
  Y(4) = X(4) + (alpha*(1.0-X(4)) - beta*X(4))*dt;
      
  % -- n1 --
  alpha = 0.014*(-44.0-X(1))/(exp(-(44.0+X(1))/2.3)-1.0);
  beta = 0.0043/exp((44.0+X(1))/34.0); 
  Y(5) = X(5) + (alpha*(1.0-X(5)) - beta*X(5))*dt;
      
      