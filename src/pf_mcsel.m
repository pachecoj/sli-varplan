function [x,Ks,y,w,w_unnorm] = pf_mcsel( ...
  num_particles, resamp_int, loc, y_all, mu_0, sig_0, sig_u, alpha  )
% PF_MCSEL - Performs sensor selection using particle filter inference and
%   particle approximation of mutual information
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
            
    % estimate mutual information
    MI = NaN(N_sensors,1);
    for k=1:N_sensors
      
      % sample Y_t
      mu_y = x_sim.^2./20;
      sig_y = sqrt( alpha * ( loc(k) - x_sim + eps).^2 );
      y_sim = mu_y + sig_y .* randn([num_particles,1]);
      
      % for each particle
      logP = NaN(num_particles,1);
      for i=1:num_particles
        this_x = x_sim(i);
        this_y = y_sim(i);
        
        % compute marginal p(x_t | y_1^{t-1} )
        tmpLogPx = -0.5*log(2*pi*sig_u*sig_u) - 0.5*(this_x - mu_sim).^2./sig_u./sig_u;
        Px = w(:,t-1)'*exp(tmpLogPx);
        
        % compute marginal p(y_t | y_1^{t-1} )
        tmpLogPy = -0.5*log(2*pi.*sig_y.*sig_y) - 0.5 * (this_y - mu_y).^2./sig_y./sig_y;
        Py = w(:,t-1)'*exp(tmpLogPy);
        
        % compute p(x_t, y_t | y_1^{t-1} )        
        tmpLogPxy = tmpLogPx - 0.5*log(2*pi.*sig_y(i).*sig_y(i)) - 0.5 * (this_y - mu_y(i)).^2./sig_y(i)./sig_y(i);
        Pxy = w(:,t-1)'*exp(tmpLogPxy);
        
        logP(i) = log(Pxy) - log(Py); % - log(Px) 
      end
      MI(k) = w(:,t-1)' * logP;
    end 
    
    % do selection
    [~,Ks(t)] = max(MI);
    y(t) = y_all(t,Ks(t));
    
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


