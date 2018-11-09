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
function U=velocity(X,kappa,kb,ks,ds,ep,mu);
  
  % reshape the position vector
  %
  N = length(X)/2;
  X  = reshape(X,N,2);
  
  % form the stokelet matrix
  %
  M = form_reg_stokes_matrix(X,ep,mu);
  e = ones(N,1);
  e(1)=0.5; e(N)=0.5;
  M = ds*M*diag([e;e]);
  
  % Evaluate the forces at the current position
  %
  [Fb,Kx] = bending_force_vec(X,kappa,kb,ds);  
  [Fs,St] = stretch_force_vec(X,ks,ds);
  F  = Fb + Fs;
  
  % compute the velocity
  %
  U = M*F(:);
  U = U(:);
  
