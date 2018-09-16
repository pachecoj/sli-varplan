function [marg, a, b] = sum_product_hmm( prior, trans, like )
% SUM_PRODUCT_HMM - Discrete sum-product on an HMM.
%
% J. Pacheco, 2017
%

  T = numel(like);

  % initialize
  Ngrid = numel(prior);
  a = zeros(Ngrid,T);
  b = zeros(Ngrid,T);

  %
  % FORWARD PASS
  %
  for t=1:T
    if t>1
      preM = a(:,t-1) + like{t-1};
      preM = exp( trans(:,:,t) + repmat(preM,[1,Ngrid,1] ) );
      a(:,t) = log( sum( preM, 1 ) );
      if any( isnan( a(:,t) ) ), keyboard; end
    else
      a(:,1) = prior;
    end
    
    % normalize message
    tmp_a = a(:,t) - max( a(:,t) );
    a(:,t) = tmp_a - log(sum(exp(tmp_a)));
  end
  
  %
  % BACKWARD PASS
  %  
  for t=(T-1):-1:1
    preM = b(:,t+1) + like{t+1};
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
    tmp_marg = a(:,t) + b(:,t) + like{t};
    tmp_marg = tmp_marg - max(tmp_marg);
    marg(:,t) = tmp_marg - log(sum(exp(tmp_marg)));
  end
  marg(:,end) = a(:,T);
  