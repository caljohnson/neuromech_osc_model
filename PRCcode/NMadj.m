function [Z2] = NMadj(Z,XLC)
% function - takes time step for adjoint problem for NM model

  global t_f t_n t_m c I dt a thresholding_on
  
  J = zeros(5,5); %store Jacobian of ODE system
  if thresholding_on == 1
    J(1,:) = (1/t_f).*[-1, 0, 0, -c*(sech(XLC(4) -2))^2, c*(sech(XLC(5) -2))^2];
  else
    J(1,:) = (1/t_f).*[-1, 0, 0, -c, c];
  end
  J(2,:) = (1/t_n).*[1, 1-3*a*XLC(2)^2, 0, 0, 0];
  J(3,:) = (1/t_n).*[-1, 0, 1-3*a*XLC(3)^2, 0, 0];
  J(4,:) = (1/t_m).*[0, 1, -1, -1, 0];
  J(5,:) = (1/t_m).*[0, -1, 1, 0, -1];
    
  dZ = -J'*Z;
  
  Z2=Z+dZ(1:5)*(-dt);

 
  
      