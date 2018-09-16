function [marg, K, a, b] = sum_product_hmm_viselect( x, prior, trans, like_all, loc, alpha, sig_x )
% SUM_PRODUCT_HMM_VISELECT - Discrete sum-product on an HMM with
%   sensor selection performed by variational information maximmization.
%
% J. Pacheco, 2017
%

  T = numel(like_all);
  K = zeros(T,1);

  % initialize
  Ngrid = numel(prior);
  a = zeros(Ngrid,T);
  b = zeros(Ngrid,T);

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
    K(t) = vi_select(x, exp(a(:,t)), loc, alpha, sig_x);    
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
  