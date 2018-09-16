function [x,Ks,y,w,w_unnorm] = pf_mcsel_lorentz( ...
  num_particles, resamp_int, F, N, deltaT, y_all, mu_0, sig_0, sig_x, sig_y  )
% PF_MCSEL_LORENTZ - Performs sensor selection using particle filter inference and
%   particle approximation of mutual information for the Lorentz95 model.
%
% J. Pacheco, 2017
%

  T = size(y_all,1);
  w = zeros([num_particles, T]);
  w_unnorm = zeros([num_particles, T]);
  x = zeros([num_particles, N, T]);
  Ks = zeros(T,1);
  N_sensors = numel( sig_y );
  y = NaN(T,N);
  odefun = @(t,z) lorentz(t,z,F,N);
  
  % sample first sensor uniformly
  Ks(1) = 1;
  y(1,:) = y_all( 1, :, Ks(1) );

  % sample initial particles
  x_sim = mu_0 + sig_0 * randn([num_particles,N]);
  D = bsxfun(@minus, x_sim, y(1,:));
  D2 = 0.5 * D.^2 ./ sig_y(Ks(1)) ./ sig_y(Ks(1));
  log_w_sim = sum(-0.5*log(2*pi.*sig_y(Ks(1)).*sig_y(Ks(1))) - D2, 2);  
  log_w_sim = log_w_sim - max(log_w_sim);
  w(:,1) = exp( log_w_sim - log(sum(exp(log_w_sim))) );
  w_unnorm(:,1) = exp( log_w_sim );
  
  % remaining particles
  for t=2:T
        
    % resample?
    if ~mod(t, resamp_int)
      cdf = cumsum( w(:,t-1) );
      [~,I] = histc(rand(num_particles,1), [0;cdf]);
      x(:,:,t-1) = x_sim( I,: );
      w_unnorm(:,t-1) = w_unnorm(I,t-1);
      w(:,t-1) = 1 / num_particles;      
    else
      x(:,:,t-1) = x_sim;
    end
       
    % solve diffeqs
    f = NaN(num_particles, N);
    for p=1:num_particles
      [~,this_f] = ode23(odefun, [0, deltaT], x(p,:,t-1));
      f(p,:) = this_f(end,:);
    end
    
    % propagate particles    
    x_sim = f + sig_x * randn([num_particles,N]);
        
    % do selection
    Ks(t) = 1;
    y(t,:) = y_all(t,:,Ks(t));
    
    % update importance weights
    D = bsxfun(@minus, x_sim, y(t,:));
    D2 = 0.5 * D.^2 ./ sig_y(Ks(t)) ./ sig_y(Ks(t));
    log_w_sim = sum(-0.5*log(2*pi.*sig_y(Ks(t)).*sig_y(Ks(t))) - D2, 2);
    log_w_sim = log_w_sim - max(log_w_sim);
    w(:,t) = exp( log_w_sim - log(sum(exp(log_w_sim))) );  
    w_unnorm(:,t) = exp( log_w_sim );
    
  end
  
  x(:,:,T) = x_sim;  
end


