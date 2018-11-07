function [Z2] = NMadj(Z,XLC)
% function - calculates time step for adjoint problem for NM model

  global t_f t_n t_m c I dt
  
  J = zeros(5,5); %store Jacobian of ODE system
  J(1,:) = (1/t_f).*[-1, 0, 0, c*(sech(XLC(4) -2))^2, -c*(sech(XLC(5) -2))^2];
  J(2,:) = (1/t_n).*[-1, 1-3*XLC(2)^2, 0, 0, 0];
  J(3,:) = (1/t_n).*[1, 0, 1-3*XLC(3)^2, 0, 0];
  J(4,:) = (1/t_m).*[0, 1, -1, -1, 0];
  J(5,:) = (1/t_m).*[0, -1, 1, 0, -1];
    
  dZ = -J'*Z;
  
  Z2=Z+dZ(1:5)*(-dt);

 
  
      