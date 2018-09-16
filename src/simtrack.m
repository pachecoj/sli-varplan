function [x, y] = simtrack( T, loc, alpha, mu_0, sig_0, sig_u)
% SIMDATA - Simulates the nonlinear state-space model in Example 1 of 
%   O. Cappe, S. Godsill, & E. Moulines, "An Overview of Existing Methods 
%     and Recent Advances in Sequential Monte Carlo". Proceedings of the 
%     IEEE, vol. 95, pp. 899-924, 2007
%
% INPUTS:
%   Ks - Sensor selections (length = #-scans)
%   loc - vector of sensor locations (length = #-sensors)
%   alpha - Measurement variance factor
%   mu_0, sig_0 - Prior mean/stdev for latent state
%   sig_u - Process noise stdev
%
% Author: J. Pacheco, 2017
%  
  x = zeros(T,1);
  N_sensors = numel(loc);
  y = zeros(T,N_sensors);
  loc = reshape(loc, [1, N_sensors]);
  
  % sample model
  for t=1:T
    
    % sample state
    if t==1
      x(1) = mu_0 + sig_0 * randn();
    else
      mu_x = x(t-1)/2 + 25*x(t-1)/(1+x(t-1)^2) + 8*cos(1.2*t);
      x(t) = mu_x + sig_u*randn();
    end
    
    % sample observations
    sig_y = sqrt( alpha * ( loc - x(t) ).^2 );
    y(t,:) = x(t)^2/20 + sig_y * randn(); 
  end
end