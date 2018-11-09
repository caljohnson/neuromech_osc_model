%
% given a vector of positions stored as X=[x;y]
%  take a time step of 
%    X(n+1) - X(n)
%    ------------- = U( X(n+1) )
%          dt
%
%  X = position of the rod
%  kappa = target curvature
%  kb = bending stiffness
%  ks = stretching stiffness
%  ds = point spacing
%  ep = regularization parameter for stokelets
%
function Xnew=BE_step3D(X,kappa,dt,kb,ks,ds,ep,mu);

  % record the number of unknowns
  %
  N = size(X,1);
    
  
  G = @(Y)( zerothis(Y,X(:),kappa,dt,kb,ks,ds,ep,mu) );
  opts = optimoptions('fsolve','SpecifyObjectiveGradient',true,'Display','none');
  Xnew = fsolve(G,X(:),opts);
  Xnew = reshape(Xnew,N,3);
  
   
  
function [Z,J]=zerothis(Xnew,X,kappa,dt,kb,ks,ds,ep,mu);
  
  Z = (Xnew(:)-X(:))/dt - velocity3D(Xnew(:),kappa,kb,ks,ds,ep,mu);
  J=velocity_jac3D(Xnew(:),kappa,kb,ks,ds,ep,mu);
  J = eye(length(X(:)))/dt - J;
  
