%
% given a vector of positions stored as X=[x;y]
%  compute the velocity on the curve
%
%  X = position of the rod
%  kappa = target curvature
%  kb = bending stiffness
%  ks = stretching stiffness
%  ds = point spacing
%  ep = regularization parameter for stokelets
%
function U=velocity3D(X,kappa,kb,ks,ds,ep,mu);
  
  % reshape the position vector
  %
  N = length(X)/3;
  X  = reshape(X,N,3);
  
  % form the stokelet matrix
  %
  M = form_reg_stokes_matrix_3D(X,ep,mu);
  e = ones(N,1);
  e(1)=0.5; e(N)=0.5;
  M = ds*M*diag([e;e;e]);
  
  % Evaluate the forces at the current position
  %  assume displacement in Z-direction is zero
  %
  [Fb,Kx] = bend_force_vec(X(:,1:2),kappa,kb,ds);  
  [Fs,St] = stretch_force_vec(X(:,1:2),ks,ds);
  F  = Fb + Fs;
  F  = [F, zeros(N,1)];
  
  % compute the velocity
  %
  U = M*F(:);
  U = U(:);
  
