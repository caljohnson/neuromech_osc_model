%
% given a vector of positions stored as X=[x;y]
%  compute the Jacibian of the velocity vector
%
%  X = position of the rod
%  kappa = target curvature
%  kb = bending stiffness
%  ks = stretching stiffness
%  ds = point spacing
%  ep = regularization parameter for stokelets
%
function J=velocity_jac3D(X,kappa,kb,ks,ds,ep,mu);
  
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
 
  % compute the force Jacobians
  %
  Jb = bend_force_jac(X(:,1:2),kappa,kb,ds);
  Js = stretch_force_jac(X(:,1:2),ks,ds);
  JF = Jb + Js;

  JF = blkdiag(JF,eye(N));
  
  % make the jacobian
  %
  J = M*JF;
  