function [K, VI] = vi_select(x, w, loc, alpha, sig_x)
% VI_SELECT - Perform sensor selection using variational information
%   maximization.
%

  N_sensors = numel( loc );    
  VI = NaN( N_sensors, 1 );
    
  x = reshape(x, [1, numel(x)]);
  w = reshape(w, [1, numel(w)]);
  
  for k=1:N_sensors
    
    % likelihood moments
    this_loc = loc(k);
    sig_y_sq = alpha * ( this_loc - x + eps).^2;
    mu_y = x.^2./20;
    
    % "matricized" quantities
    M = mu_y.*w./sig_x./sig_x;
    Lambda_k = (1/sig_x/sig_x) * w * ( sig_y_sq + mu_y.*mu_y )';    
    VI(k) = 0.5 * trace(M * (x') * x * M' / Lambda_k );
    
  end
  
  [~,K] = max(VI);
end