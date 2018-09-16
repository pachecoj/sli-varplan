function [x,Ks,y,w,w_unnorm] = pf_varinf( ...
  num_particles, resamp_int, loc, y_all, mu_0, sig_0, sig_u, alpha  )
% PF_VARINF - Performs sensor selection using particle filter inference and
%   variational information maximization on a nonlinear state-space model.
%
% J. Pacheco, 2017
%
  T = size(y_all,1);
  w = zeros([num_particles, T]);
  w_unnorm = zeros([num_particles, T]);
  x = zeros([num_particles, T]);
  Ks = zeros(T,1);
  N_sensors = numel( loc );
  y = NaN(T,1);
  
  % sample first sensor uniformly
  Ks(1) = randi( N_sensors );
  y(1) = y_all( 1, Ks(1) );

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
    
    % CHEAT: discretize p(x_t | y_1, ..., y_t )
    x_grid = linspace(-50, 50, 1000);
    mu = x_sim/2 + 25*x_sim./(1+x_sim.^2) + 8*cos(1.2*t);
    D = pdist2(x_grid', mu).^2;
    p = 1/sqrt(2*pi*sig_u*sig_u) * exp(-1/2 * D ./sig_u./sig_u);
%     dx = x_grid(2) - x_grid(1);
    p_xt = p * w(:,t-1);
    p_xt = p_xt ./ sum(p_xt);
    
    % do sensor selection    
    [Ks(t), dbg_VI] = vi_select(x_grid, p_xt, loc, alpha, 1);
    y(t) = y_all(t,Ks(t));
    
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


