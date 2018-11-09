%
% compute the curvature of a filament in the plane
%   curvature is computed by using that
%     k = N.dT/ds 
%   discretely this is
%     k = (Np+Nm)/2 * (Tp-Tm)/ds
%   Where Tp and Tm are the forward and backward differences
%   the above formula simplifes to 
%     k = (Tm x Tp).Z/ds
%   Use a coordinate system on the filment so that
%     T x N = Z
%
%   This code assumes that length of each segment is (approximately)
%   ds. Thus the tangent vector is computed as the difference of
%   adjacent points on the filament.
%
%   **Note that curvature is only computed at the interior points. We
%   are treating the curve as piecewise linear, and computing the
%   change in tangent angle at interior points. 
%
function K = curvature(X,ds);
  
  % record the number of points
  %
  N = size(X,1);
  
  % compute the differences of the point location
  %  note that D(i) is the forward difference for point i
  %
  D = X(2:N,:) - X(1:N-1,:);
  Dp = [D; [0 0]];
  Dm = [[0 0]; D];

  % compute the curvature
  %
  K = zeros(N,1);
  J = 2:N-1; 
  K(J) =  (Dp(J,2).*Dm(J,1) - Dp(J,1).*Dm(J,2))/ds^3;
  
  