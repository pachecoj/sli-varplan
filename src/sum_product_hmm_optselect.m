function [marg, K, a, b] = sum_product_hmm_optselect( x_grid, prior, trans, like_all, loc, alpha, sig_x )
% SUM_PRODUCT_HMM_OPTSELECT - Discrete sum-product on an HMM with
%   optimal sensor selection performed by numerical integration.
%
% J. Pacheco, 2017
%

  T = numel(like_all);
  K = zeros(T,1);
  N_sensors = numel(loc);
  small_logval = -1000;

  % initialize
  Ngrid = numel(prior);
  a = zeros(Ngrid,T);
  b = zeros(Ngrid,T);
  
  % discrete grid on XY
  y_grid = linspace(-200, 200, 2*Ngrid);
  [X,Y] = meshgrid( x_grid, y_grid );
  xy_grid = [ X(:) Y(:) ];
  dx = x_grid(2) - x_grid(1);
  dy = y_grid(2) - y_grid(1);
%   clear X Y

  %
  % FORWARD PASS
  %
  for t=1:T
    if t>1
      preM = a(:,t-1) + like_all{t-1}(:,K(t-1));
      preM = exp( trans(:,:,t) + repmat(preM,[1,Ngrid,1] ) );
      a(:,t) = log( sum( preM, 1 ) );
      if any( isnan( a(:,t) ) ), keyboard; end
    else
      a(:,1) = prior;
    end
        
    % normalize message
    tmp_a = a(:,t) - max( a(:,t) );
    a(:,t) = tmp_a - log(sum(exp(tmp_a)));
    
    % do sensor selection
    log_prior = a(:,t) - log(dx);
    Hx =  - exp(log_prior)' * log_prior * dx;
    MI = NaN(N_sensors,1);    
    for l=1:N_sensors
      
      % discretize likelihood
      sig_y = sqrt( alpha * ( loc(l) - xy_grid(:,1) + eps).^2 );
      mu_y = xy_grid(:,1).^2./20;
      ll = -0.5 * log(2*pi*sig_y.*sig_y) - 0.5 * (xy_grid(:,2) - mu_y).^2 ./ sig_y ./ sig_y;
      ll = reshape( ll, [2*Ngrid, Ngrid] )'; % X-by-Y
      ll = bsxfun(@minus, ll, max( ll, [], 2 ));
%       ll = bsxfun(@minus, ll, log(sum(exp(ll),2))) - log(dy);
      
      % compute marginal log p(y) and normalize
      log_joint = bsxfun(@plus, ll, log_prior - max(log_prior));
      log_joint = log_joint - max(log_joint(:));
      log_yMarg = log(sum(exp(log_joint), 1));
      log_yMarg = log_yMarg - max(log_yMarg);
      log_yMarg( log_yMarg(:) < small_logval ) = small_logval;
%       log_yMarg = log_yMarg - log(sum(exp(log_yMarg))) - log(dy);
      
      % compute posterior log p(x|y) and normalize
      log_post = bsxfun(@minus, log_joint, log_yMarg);
      log_post = bsxfun(@minus, log_post, max(log_post,[],1));
      log_post( log_post(:) < small_logval ) = small_logval;
      log_post = bsxfun(@minus, log_post, log(sum(exp(log_post),1))) - log(dx);
            
%       p_joint = exp(log_joint) ./ sum(exp(log_joint(:))) ./ dx ./ dy;
%       p_post = bsxfun(@rdivide, p_joint, sum(p_joint,1)*dx);      
%       log_p_post = log(p_post);
%       log_p_post( log_p_post(:) < small_logval ) = small_logval;
      
      % compute mutual information
      p_joint = exp(log_joint) ./ sum(exp(log_joint(:))) ./ dx ./ dy;
      H_post = - p_joint(:)' * log_post(:) * dx * dy;
      MI(l) = Hx - H_post;   
      
%      % DEBUG:      
%      figure('InvertHardcopy','off','Color',[1 1 1]);
%      set(gcf,'Position',[200 350 1550 400]);
%      colormap gray
%      subplot(1,3,1);
%      contour(X,Y,exp(ll)');     
%      hold on;
%      plot([loc(l) loc(l)], [y_grid(1) y_grid(end)], '-r');
%      title(sprintf('Likelihood P(Y|X)',loc(l)));
%      xlabel('Position (X)');
%      ylabel('Observation (Y)');
%      legend('PDF','Sensor Location');
%      subplot(1,3,2);
%      contour(X,Y,exp(log_joint)');
%      hold on;
%      plot([loc(l) loc(l)], [y_grid(1) y_grid(end)], '-r');
%      title('Joint P(X,Y)');
%      xlabel('Position (X)');
%      subplot(1,3,3);
%      contour(X,Y,exp(log_post)');
%      hold on;
%      plot([loc(l) loc(l)], [y_grid(1) y_grid(end)], '-r');
%      title('Posterior P(X|Y)');
%      xlabel('Position (X)');

    end
    [~,K(t)] = max(MI);
  end
  
  %
  % BACKWARD PASS
  %  
  for t=(T-1):-1:1
    preM = b(:,t+1) + like_all{t+1}(:,K(t+1));
    preM = exp( trans(:,:,t+1) + repmat(preM', [Ngrid,1,1]) );
    b(:,t) = log( sum( preM, 2 ) );
    
    % normalize
    tmp_b = b(:,t) - max( b(:,t) );
    b(:,t) = tmp_b - log(sum(exp(tmp_b)));
  end  
  
  %
  % COMPUTE  MARGINALS
  %
  marg = zeros(Ngrid,T);    
  for t=1:(T-1)
    tmp_marg = a(:,t) + b(:,t) + like_all{t}(:,K(t));
    tmp_marg = tmp_marg - max(tmp_marg);
    marg(:,t) = tmp_marg - log(sum(exp(tmp_marg)));
  end
  marg(:,end) = a(:,T);
  