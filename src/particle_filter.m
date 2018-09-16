function [x,w,w_unnorm] = particle_filter( ...
  num_particles, resamp_int, Ks, loc, y, mu_0, sig_0, sig_u, alpha  )
% PARTICLE_FILTER - Runs standard particle filter on nonlinear state-space
%   model of Example 1 in Cappe et. al.
%
% J. Pacheco, 2017
%
  T = numel(y);
  w = zeros([num_particles, T]);
  w_unnorm = zeros([num_particles, T]);
  x = zeros([num_particles, T]);

  % sample initial particles
  x_sim = mu_0 + sig_0 * randn([num_particles,1]);
  mu_y = x_sim.^2./20;
  sig_y = sqrt( alpha * ( loc(Ks(1)) - x_sim + eps).^2 );
  log_w_sim = -0.5*log(2*pi.*sig_y.*sig_y) - 0.5*(y(1) - mu_y).^2./sig_y./sig_y;  
  log_w_sim = log_w_sim - max(log_w_sim);
  w(:,1) = exp( log_w_sim - log(sum(exp(log_w_sim))) );
  w_unnorm(:,1) = exp( log_w_sim );
  
  % remaining particles
  for t=2:T
    
    % resample?
    if ~mod(t, resamp_int)
      cdf = cumsum( w(:,t-1) );
      [~,I] = histc(rand(num_particles,1), [0;cdf]);
      x(:,t-1) = x_sim( I );
      w_unnorm(:,t-1) = w_unnorm(I,t-1);
      w(:,t-1) = 1 / num_particles;      
    else
      x(:,t-1) = x_sim;
    end
    
    % propagate
    mu_sim = x(:,t-1)/2 + 25*x(:,t-1)./(1+x(:,t-1).^2) + 8*cos(1.2*t);
    x_sim = mu_sim + sig_u * randn([num_particles,1]);
    
    % update importance weights
    mu_y = x_sim.^2./20;
    sig_y = sqrt( alpha * ( loc(Ks(t)) - x_sim + eps).^2 );
    log_w_sim = log(w(:,t-1)) -0.5*log(2*pi.*sig_y.*sig_y) ...
      - 0.5*(y(t) - mu_y).^2./sig_y./sig_y; 
    log_w_sim = log_w_sim - max(log_w_sim);
    w(:,t) = exp( log_w_sim - log(sum(exp(log_w_sim))) );  
    w_unnorm(:,t) = exp( log_w_sim );
    
  end
  
  x(:,T) = x_sim;  
end