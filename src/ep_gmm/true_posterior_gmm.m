function [ true_posterior_pdf, logZ, true_posterior_gauss ] = true_posterior_gmm( ...
  gmm, ...    % GMM parameters
  y, ...      % Data vector
  x_vals, ... % Discrete grid
  k_vec ...   % Sensor selection
  )
% TRUE_POSTERIOR_GMM - Computes the true posterior PDF for observations
%   drawn iid from a two-component GMM with specified parameters.
%
% RETURNS:
%   true_posterior_pdf - Discretized PDF computed on grid values 'x_vals'
%   logZ - Log-partition (log-likelihood)
%   true_posterior_gauss - Best Gaussian approximation to posterior
%
% Jason L. Pacheco
% 12/15/11
%
  N = numel(y);
  dx = x_vals(2) - x_vals(1);
  
  % true posterior
  px = normpdf( x_vals, gmm.m0_x, sqrt(gmm.V0_x) )';
  p_yGx = zeros( numel(x_vals), N );
  for n = 1:N
    k = k_vec(n);
    y_i = y(n);
    p_yGx(:,n) = dense_mixture( gmm, x_vals, y_i, k ); 
  end
  true_posterior_pdf = ( px .* prod( p_yGx, 2 ) )';
  Z = sum( true_posterior_pdf ) * dx;
  true_posterior_pdf = true_posterior_pdf ./ Z;
  logZ = log(Z);    
  
  % best Gaussian approx
  Ex = x_vals * true_posterior_pdf' * dx;
  Exx = x_vals.^2 * true_posterior_pdf' * dx;
  Var = Exx - Ex.^2;
  true_posterior_gauss = normpdf(x_vals, Ex, sqrt(Var));
end

