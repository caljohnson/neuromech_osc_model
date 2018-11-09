function [ driving_moment ] = moment_step( X,c,no_mods, S)
%MOMENT_STEP Creates the driving moment to the beam
%   

%set vector of Av inputs
AV = zeros(no_mods,1);
AD = zeros(no_mods,1);
for jj=0:no_mods-1
  AV(jj+1) = X(3+jj*4);
  AD(jj+1) = X(4+jj*4);
end

%spreads output of each moment to section of beam
driving_moment = c*(tanh(S*AV-2) - tanh(S*AD-2));

end

