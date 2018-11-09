%
% Compute the bending force per unit length along a 2D inextensible
% rod.
%
% X       array with positions of the rod 
% kappa   resting curvature where there is no force
% kb      bending stiffness
% ds      spacing between the points
%
% The force is computed from taking the variation of the discrete
% energy
%    E = kb/2 * sum( (Kx-K0)^2, j=2,N-1 ) * ds
% where Kx is the current curvature at the interior nodes.
%
%   curvature is computed by using that
%     Kx = N.dT/ds 
%   discretely this is
%     Kx = (Np+Nm)/2 * (Tp-Tm)/ds
%   Where Tp and Tm are the forward and backward differences
%   the above formula simplifes to 
%     Kx = (Tm x Tp).Z/ds
%   Use a coordinate system on the filment so that
%     T x N = Z
%   It is assumed that the points are spaced by ds so that the tangent
%   vector is just Tp = (X(j+1)-X(j))/ds
%
function [F,Kx] = bend_force(X,kappa,kb,ds);
  
  % record the number of points
  %
  N = size(X,1);
  
  % initialize the forces
  %
  F = zeros(N,2);
  
  
  % compute the differences of the point location
  %  note that D(i) is the forward difference for point i
  %
  D = X(2:N,:) - X(1:N-1,:);
  Dp = [D; [0 0]];
  Dm = [[0 0]; D];
  
  % compute the energy density at each point
  %
  Kx = zeros(N,1);
  W = zeros(N,1);
  K = 2:N-1;
  Kx(K) =  (Dp(K,2).*Dm(K,1) - Dp(K,1).*Dm(K,2))/ds^3;
  W(K) = Kx(K) - kappa(K);
  
  % update the bending forces
  %
  J = 2:N-1;
  WW = repmat(W,1,2);
  F(J-1,:) = F(J-1,:) - WW(J,:).*[ -Dp(J,2)        ,  Dp(J,1)        ];
  F(J  ,:) = F(J  ,:) - WW(J,:).*[  Dm(J,2)+Dp(J,2), -Dp(J,1)-Dm(J,1)];
  F(J+1,:) = F(J+1,:) - WW(J,:).*[ -Dm(J,2)        ,  Dm(J,1)        ];
 
  
% $$$   
% $$$   
% $$$   % loop over the interior points and update forces
% $$$   %
% $$$   G = zeros(N,2);
% $$$   for j=2:N-1
% $$$    G(j-1,:) = G(j-1,:) - W(j)*[  -Dp(j,2)        ,  Dp(j,1)        ];
% $$$    G(j  ,:) = G(j  ,:) - W(j)*[   Dm(j,2)+Dp(j,2), -Dp(j,1)-Dm(j,1)];
% $$$    G(j+1,:) = G(j+1,:) - W(j)*[  -Dm(j,2)        ,  Dm(j,1)        ];
% $$$   end
% $$$   fprintf('max diff of F and G %g \n',max( abs(F(:)-G(:))));  
% $$$   
  
  % rescale the forces
  %
  F = kb * F/ds^3;
