function [ep_state, skipped_any] = ep_gmm(...
  gmm, ...  % Gaussian mixture model parameters
  opt, ...  % EP options
  x, ...    % true latent variable (for computing error)
  y, ...    % data vector
  k_vec, ... % Sensor selection
  x_grid, ... % X grid points (for plotting)  
  fnplot, ... % Plotting function
  h_fig ... % figure handle
  )
% EP_GMM - implements Expectation Propagation for a model with i.i.d.
% Gaussian Mixture factors.  The likelihood factors are of the form:
%
% p(y|x) \propto w * N(y;Ax+a,V_0) + (1-w) N(y;Bx+b,V_1)
%
% RETURNS:
%   ep_state  - State of EP including posterior moments/messages/augmented
%               dist. moments.
%   err - error vector for each iteration
%
% Jason L. Pacheco
% 9/3/10
%

  N = numel(y);
  dx = x_grid(2) - x_grid(1);

  %
  % CREATE FIGURE HANDLES
  %
  if opt.plot_flag && ~isempty(h_fig)
    figure(h_fig);
    h_plot = subplot(2,1,1);
    h_err = subplot(2,1,2);     
    subplot(h_plot), cla;
    subplot(h_err), cla;    
  end
  
%   %
%   % CALCULATE TRUE POSTERIOR
%   %
%   x_vals = linspace(-100, 100, 10000 );
%   dx = x_vals(2) - x_vals(1);
%   true_posterior_pdf = true_posterior_gmm( gmm, y, x_vals );
%   m_x_true = sum( x_vals .* true_posterior_pdf ).*dx;
%   v_x_true = sum( (x_vals - m_x_true ).^2 .* true_posterior_pdf ) .* dx;
%   best_fit_pdf = normpdf( x_vals, m_x_true, sqrt(v_x_true) );

  % Initialize terms
  eta = zeros(1,N);
  lam = zeros(1,N);
  m_x = gmm.m0_x;
  v_inv_x = 1/gmm.V0_x;
  v_x = 1/v_inv_x;
  msg_pdf = zeros(numel(x_grid),N);
  V_xa = zeros([N,1]); V_xb = zeros([N,1]); 
  m_xa = zeros([N,1]); m_xb = zeros([N,1]);

  %
  % EP - MAIN LOOP
  %
  iters = 0; done = false; 
  skipped_any = false;
  while ~done

    iters = iters+1;
    lam_old = lam;
    eta_old = eta;
    
    % iterate over factors
    skipped_this = false;
    for i=opt.sched %1:N      

      % update cavity distribution
      [ v_inv_divi, m_divi ] = ...
        cavity_gauss( v_inv_x, m_x, lam(i), eta(i) );

      % get posterior moments
      try     
        p_aug = dense_augmented_gmm( ...
          gmm, x_grid, y(i), m_divi, 1/v_inv_divi, k_vec(i) );          
        [ m_x_aug, ~, v_inv_x_aug ] = dense_moments( x_grid, p_aug );
      catch ex
        if (strcmp(ex.identifier, 'UpdateAugmentedGMM:MixWts'))
          skipped_this = true;
          continue;
        else
          rethrow(ex)
        end
      end
      
      % set stepsize
      switch(opt.damp_type)
        case 0
          this_damp_fact = 1.0;
        case 1
          this_damp_fact = opt.min_damp_fact;
      end      

      % update posterior moments
      [v_inv_x_new, m_x_new] = ...
        damp_posterior( v_inv_x_aug, m_x_aug, v_inv_x, m_x, this_damp_fact );
      if ( v_inv_x_new > 0.0 )
        v_inv_x = v_inv_x_new;
        m_x = m_x_new;
      else
        skipped_this = true;
      end       
      
      % update msg parameters
      [lam(i), eta(i)] = ...
        msg_update_gauss( v_inv_x, m_x, v_inv_divi, m_divi );
      
      % update msg PDF
      msg_pdf(:,i) = exp(-0.5 * x_grid.^2 * lam(i) + x_grid .* eta(i));
      msg_pdf(:,i) = msg_pdf(:,i) ./ sum( msg_pdf(:,i) ) ./ dx;
      
    end
    skipped_any = skipped_this | skipped_any;

    %
    % TEST CONVERGENCE
    %
    if iters == opt.max_iters
      done = true;
    elseif (iters > 1) && (~skipped_this)
      done = test_conv( lam, eta, lam_old, eta_old, opt.conv_thresh );      
    end

    %
    % CALCULATE ERROR
    % L1 distance between approx. and true posterior
    %
    y_hat_pdf = 1/sqrt(2*pi*v_x) ...
      .* exp( -0.5 .* (x_grid - repmat(m_x,size(x_grid)) ).^2.*v_inv_x );  
    if sum( y_hat_pdf ) > 0      
      y_hat_pdf = y_hat_pdf / sum( y_hat_pdf ) / dx;
    else
      y_hat_pdf = 1/numel(y_hat_pdf(:)) + y_hat_pdf;
    end
%     err(end+1) = l1_diff_symmetric( y_hat_pdf, true_posterior_pdf, x_vals );
%     if isnan(err(end)) , keyboard; end;    
%     fprintf('Iteration #%d: %E\n', iters, err(end) );

    %
    % PLOT ESTIMATES
    %
    if opt.plot_flag && ~isempty(h_fig)
%       plot_ep_iter( h_fig, true_posterior_pdf, y_hat_pdf, ...
%         best_fit_pdf, msg_pdf, err, y, x_vals );
      fnplot(h_fig, y_hat_pdf, msg_pdf, y, x_grid)
    end

  end
  
  %
  % SAVE EP STATE
  %
  ep_state = struct('v_inv_x', v_inv_x, 'm_x', m_x, 'lam', lam, 'eta', eta, ...
    'V_xa', V_xa, 'V_xb', V_xb, 'm_xa', m_xa, 'm_xb', m_xb, 'iters', iters);

