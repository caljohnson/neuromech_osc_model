function [ Y ] = NM_step3d( X, K, Pa, dt, no_mods)
%NM_STEP3D Takes a forward euler step for each neuromechanical oscillator
%unit
%   

%module parameters
a = 1;
t_n = 0.1;
t_m = 1;
I = 0;

%interpolate curvature input to each module
K = Pa*K;
Y = zeros(4*no_mods,1);
for jj=0:no_mods-1
 % -- V1 --
  Y(1+jj*4) = X(1+jj*4) + (X(1+jj*4)-a*X(1+jj*4)^3 + I - K(1+jj))*dt/t_n;
  
  % -- V2 --
  Y(2+jj*4) = X(2+jj*4) + (X(2+jj*4)-a*X(2+jj*4)^3 + I + K(1+jj))*dt/t_n;
  
  % -- A1 --
  Y(3+jj*4) = X(3+jj*4) + (-X(3+jj*4) + X(1+jj*4) - X(2+jj*4))*dt/t_m;
  
  % -- A2 --
  Y(4+jj*4) = X(4+jj*4) + (-X(4+jj*4) + X(2+jj*4) - X(1+jj*4))*dt/t_m;
end

end

