function [Y] = NMchain_step(X,step)
% function - takes time step for NM chain model

  global t_f t_n t_m a_v c I dt dim A gamma W
  
%   % -- K --
%   M = tanh(X(3*dim+1:4*dim)-2) - tanh(X(4*dim+1:5*dim)-2);
%   %Forward Euler scheme
% %   dK = (gamma*speye(dim) + A*A')\(-A*A'*(X(1:dim) + c.*M)./t_f);
% %   Y(1:dim) = X(1:dim) + dK*dt;  
%     
%   %Backward Euler scheme
%   lhs_matrix = gamma/dt*speye(dim) + (1/dt + 1/t_f)*A*A';
%   rhs_stuff = (gamma*speye(dim) + A*A')*X(1:dim)/dt - c*A*A'*M/t_f;
%   Y(1:dim) = lhs_matrix\rhs_stuff;
% 
%   % -- V1 --
%   dV1 = X(dim+1:2*dim) - a_v.*X(dim+1:2*dim).^3 + I - W*X(1:dim);
%   Y(dim+1:2*dim) = X(dim+1:2*dim) + dV1*dt/t_n;
%   
%   % -- V2 --
%   dV2 = X(2*dim+1:3*dim) - a_v.*X(2*dim+1:3*dim).^3 + I + W*X(1:dim);
%   Y(2*dim+1:3*dim) = X(2*dim+1:3*dim) + dV2*dt/t_n;
%   
%   % -- A1 --
%   dA1 = -X(3*dim+1:4*dim) + X(dim+1:2*dim) - X(2*dim+1:3*dim);
%   Y(3*dim+1:4*dim) = X(3*dim+1:4*dim) + dA1*dt/t_m;
%   
%   % -- A2 --
%   dA2 = -X(4*dim+1:5*dim) - X(dim+1:2*dim) + X(2*dim+1:3*dim);
%   Y(4*dim+1:5*dim) = X(4*dim+1:5*dim) + dA2*dt/t_m;
%       

RHS_matrix = (gamma*speye(dim)+A*A')\(-A*A'/t_f);

K_dot = @(t,K, AV, AD) RHS_matrix*(K + c.*(tanh(AV-2)-tanh(AD-2)));
volt_V_dot = @(t,Kappa,volt_V) (1/t_n).*(volt_V - a_v.*volt_V.^3 + I + W*Kappa);
volt_D_dot = @(t,Kappa,volt_D) (1/t_n).*(volt_D - a_v.*volt_D.^3 + I - W*Kappa);
AV_dot = @(t,volt_V, volt_D, A) (1/t_m).*(volt_V - volt_D - A);
AD_dot = @(t,volt_V, volt_D, A) (1/t_m).*(volt_D - volt_V - A);

ode_rhss = @(t,X) [K_dot(t, X(1:dim), X(3*dim+1:4*dim), X(4*dim+1:5*dim));...
    volt_V_dot(t, X(1:dim), X(dim+1:2*dim)); ...
    volt_D_dot(t, X(1:dim), X(2*dim+1:3*dim)); ...
    AV_dot(t, X(dim+1:2*dim), X(2*dim+1:3*dim), X(3*dim+1:4*dim)); ...
    AV_dot(t, X(dim+1:2*dim), X(2*dim+1:3*dim), X(4*dim+1:5*dim)); ];
           
% options = odeset('RelTol',1e-8,'AbsTol',1e-10,  'MaxStep', max_step);
[t,y] = ode23(ode_rhss,[0,step], X);
Y = y(end,:);

      