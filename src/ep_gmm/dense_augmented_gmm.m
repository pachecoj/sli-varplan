function [pdf, Z] = dense_augmented_gmm( ...
  gmm, ...
  x_vals, ...
  y, ...
  m, ...
  v, ...
  k ...
  )
% DENSE_AUGMENTED_GMM - Multiply Gaussian prior with mixture likelihood.

  dx = x_vals(2) - x_vals(1);
  prior_pdf = normpdf( x_vals, m, sqrt(v) );
  
  % 1D posterior
  if numel(y)==1
    like_pdf = dense_mixture( gmm, x_vals, y, k );
    pdf = prior_pdf .* reshape(like_pdf,size(x_vals));
    Z = sum( pdf ) * dx;
    pdf = pdf / Z ;
    
  % 2D joint
  else
    dy = y(2) - y(1);
    like_pdf = dense_mixture_2D( gmm, x_vals, x_vals, k );
    log_pdf = bsxfun(@plus, log( like_pdf ), log( prior_pdf ));
    log_pdf = log_pdf - max(log_pdf(:));
    pdf = exp(log_pdf) ./ sum(exp(log_pdf(:))) ./ dx ./ dy;
    Z = 1;
  end
  
end