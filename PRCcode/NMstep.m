function [Y] = NMstep(X)
% function - takes time step for NM model

  global t_f t_n t_m c I dt a

  % -- K --
  Y(1) = X(1) + (-X(1) + c*(tanh(X(4)-2) - tanh(X(5) - 2)))*dt/t_f;  

  % -- V1 --
  Y(2) = X(2) + (X(2)-a*X(2)^3 + I - X(1))*dt/t_n;
  
  % -- V2 --
  Y(3) = X(3) + (X(3)-a*X(3)^3 + I + X(1))*dt/t_n;
  
  % -- A1 --
  Y(4) = X(4) + (-X(4) + X(2) - X(3))*dt/t_m;
  
  % -- A2 --
  Y(5) = X(5) + (-X(5) + X(3) - X(2))*dt/t_m;
      
      