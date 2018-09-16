function [k_opt, err, H_cond_var] = vi_plan_gmm( model, x_vals, ep_state, h_fig )
% VI_PLAN_GMM - Planning for GMM based on variational lower bound on MI.
%
% J. Pacheco, 2017

  MIN_LOG_VAL = -1e6;
  Nsensors = numel( model.loc );
  dx = x_vals(2) - x_vals(1);
    
  % best moments
  H_cond_best = Inf;  
  p_joint_best = [];
  m_best = [];
  V_best = [];
  
  H_cond_var = zeros(Nsensors,1);    
  err = NaN(Nsensors,1);
  for k=1:Nsensors
    
    % compute joint
    p_joint = dense_augmented_gmm( ...
      model, x_vals, x_vals, ep_state.m_x, 1/ep_state.v_inv_x, k );
    
    % moment-matching
    [ m_xy, V_xy ] = dense_moments( x_vals, p_joint );
    V_y = V_xy(2,2);
    m_y = m_xy(2);
    
    % compute Gaussian conditional entropy
    H_q_xy = 1/2 * log( (2*pi*exp(1))^2 * det( V_xy ) );
    H_q_y = 1/2 * log( 2*pi*exp(1)*V_y );
    H_cond_var(k) = H_q_xy - H_q_y;
    
    % compute numerical entropy
    pmarg_y = sum( p_joint, 2 ) * dx;
    logp_marg_y = log( pmarg_y );
    log_p_joint = log( p_joint );
    log_p_joint( isinf(log_p_joint) ) = MIN_LOG_VAL;
    Hxy = - p_joint(:)' * log_p_joint(:) * dx * dx;
    Hy = - pmarg_y' * logp_marg_y * dx;
    H_cond_opt = Hxy - Hy;
    err(k) = H_cond_var(k) - H_cond_opt;
      
    % best?
    if H_cond_var(k) < H_cond_best
      H_cond_best = H_cond_var(k);
      p_joint_best = p_joint;
      m_best = m_xy;
      V_best = V_xy;
      k_opt = k;
    end    
  end
  
  % plot stuff
  if ~isempty( h_fig )
    
    % compute variational distribution
    [X,Y] = meshgrid(x_vals);
    XY = [ X(:) Y(:) ];
    q_best = mvnpdf( XY, m_best', V_best );
    
    % plot
    figure(h_fig);
    subplot(2,1,2);
    cla;
    set(gca,'FontSize',14);
    set(gca,'YScale','linear');
    hold on;
    contour(X,Y,p_joint_best);    
    contour(X,Y,reshape(q_best,size(X)));
    for k=1:Nsensors
      text(model.loc(k), 0, sprintf('%d',k));
    end
    plot(model.loc(k_opt), 0, 'or', 'MarkerSize', 15)
    xlabel('X');  ylabel('Y_t');
    legend('Joint p(X,Y)','Moment-Matched q(X,Y)');
    title(sprintf('Sensor %d Joint PDF',k_opt));
    xlim([-25 25]);
    ylim([-25 25]);
  end
  
end