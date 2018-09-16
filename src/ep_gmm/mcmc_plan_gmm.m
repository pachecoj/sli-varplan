function [k_opt, err, y_samp ] = mcmc_plan_gmm(Nsamp, model, x_grid, x_samp, true_posterior_pdf, h_fig)
% MCMC_PLAN_GMM - Sensor selection based on MCMC samples.
%
% Estimates mutual information in the form:
%
% I(X;Y_{n+1} | y_1^n) = H(Y_{n+1} | y_1^n) - H(Y_{n+1} | X)
%
% J. Pacheco, 2017
%
  Nsensors = numel( model.loc );
  MI = zeros(Nsensors,1);  
  MIN_LOG_VAL = -1e6;
  dx = x_grid(2) - x_grid(1);
  
  % best estimates
  MI_best = -Inf;  
  p_joint_best = [];
  k_opt = NaN;  
  err = NaN(Nsensors,1);
  
  % loop over sensors    
  y_samp = cell(Nsensors,1);
  for k=1:Nsensors    
    fprintf('\tSensor %d of %d\n', k, Nsensors);
    mu_k = model.loc(k);
        
    % draw observation samples
    y_samp{k} = zeros(Nsamp,Nsamp);
    for i=1:Nsamp
      this_x = x_samp(i);
            
      % sample assignments
      pi_vec = rand(Nsamp,1);
      I_A = find(pi_vec <= model.w);
      I_B = find(pi_vec > model.w);
      
      % sample measurements
      dist = abs( mu_k - this_x ); 
      var_tgt = model.c * dist + model.V_1;
      y_samp{k}(i,I_A) = (model.A*this_x + model.a) + sqrt(model.V_0) * randn(1,numel(I_A));
      y_samp{k}(i,I_B) = (model.B*this_x + model.b) + sqrt(var_tgt) * randn(1,numel(I_B));
    end
    
    % conditional entropy H(Y|X)
    H_cond = 0;
    for i=1:Nsamp      
      this_x = x_samp(i);
      var_tgt = model.c * abs( mu_k - this_x ) + model.V_1;
      for j=1:Nsamp        
        this_y = y_samp{k}(i,j);
        p0 = sqrt(2*pi*model.V_0) * exp( - 0.5 * (this_y - model.A*this_x + model.a)^2 / model.V_0 );        
        p1 = sqrt(2*pi*var_tgt) * exp( -0.5 * (this_y - model.B*this_x + model.b)^2 / var_tgt );
        H_cond = H_cond - log( model.w * p0 + (1-model.w) * p1 );
      end
    end
    H_cond = (1/Nsamp/Nsamp) * H_cond;
    
    % marginal entropy H(Y)
    H_marg = 0;
    for j=1:numel(y_samp{k}(:))      
      p_hat = 0;
      for i=1:Nsamp
        this_y = y_samp{k}(j);
        this_x = x_samp(i);
        var_tgt = model.c * abs( mu_k - this_x ) + model.V_1;
        p0 = sqrt(2*pi*model.V_0) * exp( - 0.5 * (this_y - model.A*this_x + model.a)^2 / model.V_0 );        
        p1 = sqrt(2*pi*var_tgt) * exp( -0.5 * (this_y - model.B*this_x + model.b)^2 / var_tgt );
        p_hat = p_hat + model.w * p0 + (1-model.w) * p1;
      end
      log_phat = log( 1/Nsamp * p_hat );
      H_marg = H_marg - log_phat;
    end        
    H_marg = 1/numel(y_samp{k}(:)) * H_marg;
    MI(k) = H_marg - H_cond;
    
    % compute numerical conditional
    y_grid = x_grid;
    DA = pdist2( model.A.*x_grid' + model.a, y_grid');
    DB = pdist2( model.B.*x_grid' + model.b, y_grid');
    var_tgt = model.c * abs( mu_k - x_grid ) + model.V_1;
    p0 = sqrt(2*pi*model.V_0) * exp( - 0.5 * DA.^2 / model.V_0 );        
    p1_exp = bsxfun(@times, -0.5 * DB.^2, 1./var_tgt');
    p1 = bsxfun(@times, sqrt(2*pi.*var_tgt'), exp(  p1_exp ) );
    p_cond = model.w * p0 + (1-model.w) * p1;
    p_cond = bsxfun(@times, p_cond, 1./sum(p_cond,2)./dx);
    
    % compute numerical joint p_k(X,Y_{n+1} | y_1^n )    
    log_p_joint = bsxfun(@plus, log(true_posterior_pdf), log(p_cond));    
    log_p_joint = log_p_joint - max( log_p_joint(:) );
    log_p_joint( isinf(log_p_joint) ) = MIN_LOG_VAL;
    log_p_joint = log_p_joint - log( sum( exp( log_p_joint(:) ) ) ) - log(dx) - log(dx);
    p_joint = exp( log_p_joint );
    Hxy = - p_joint(:)' * log_p_joint(:) * dx * dx;
    
    % compute marginals    
    pmarg_y = sum( p_joint, 1 ) * dx;
    logp_marg_y = log( pmarg_y ) - max( log( pmarg_y ) );    
    logp_marg_y = logp_marg_y - log( sum( exp( logp_marg_y ) ) ) - log(dx);        
    Hy = - exp(logp_marg_y(:))' * logp_marg_y(:) * dx;    
    pmarg_x = sum( p_joint, 2 ) * dx;
    logp_marg_x = log( pmarg_x ) - max( log( pmarg_x ) );
    logp_marg_x = logp_marg_x - log( sum( exp( logp_marg_x ) ) ) - log(dx);    
    Hx = - exp(logp_marg_x(:))' * logp_marg_x(:) * dx;
    MI_opt = Hx + Hy - Hxy;
    err(k) = MI(k) - MI_opt;
    
    % best?
    if MI(k) > MI_best
      MI_best = MI(k);
      p_joint_best = p_joint;
      k_opt = k;
    end
  end  
  
  % plot stuff
  if ~isempty( h_fig )  
    figure(h_fig);
    subplot(2,1,2);
    cla;
    set(gca,'FontSize',14);
    set(gca,'YScale','linear');
    hold on;
    [X,Y] = meshgrid( x_grid );
    contour(X,Y,p_joint_best);    
    Xs = repmat( x_samp', [1, Nsamp] );    
    scatter(Xs(:), y_samp{k_opt}(:), '.k');
    plot(model.loc(k_opt), 0, 'or', 'MarkerSize', 15)
    for k=1:Nsensors
      text(model.loc(k), 0, sprintf('%d',k));
    end            
    xlabel('X');  ylabel('Y_t');
    legend('Joint p(X,Y_{n+1} | y_1^n)','MCMC Samples','Selected Sensor');
    title(sprintf('Sensor %d Joint PDF',k_opt));
    xlim([-25 25]);
    ylim([-25 25]);
  end
  
end