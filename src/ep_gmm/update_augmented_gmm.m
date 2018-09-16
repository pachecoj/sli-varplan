function [ V_aug, m_aug ] = update_augmented_gmm( ...
  gmm, ...        % Gauss. Mixture Model params
  v_inv_divi, ... % Cavity inverse variance
  m_divi, ...     % Cavity mean
  y_i ...         % Data for this factor
  )
% UPDATE_AUGMENTED_GMM - Updates the augmented distribution for an iid
%   two-component Gaussian mixture model.  Given mean parameters of the
%   cavity distribution, GMM parameters, and data point for this factor,
%   returns the mean, variance, and mixture weights for the GMM augmented
%   posterior.
%
% RETURNS:
%   V_aug - Variance of augmented posterior
%   m_aug - Mean of augmented posterior
%
% THROWS:
%   UpdatedAugmentedGMM:MixWts - Indicates cavity variance is sufficiently
%     negative to result in unnormalizeable mixture weights.
%
% Jason L. Pacheco
% 12/15/11
%

  v_divi = 1/v_inv_divi;

  % incorporate measurement
  V_xa = 1/( v_inv_divi + gmm.A'*1/gmm.V_0*gmm.A );
  V_xb = 1/( v_inv_divi + gmm.B'*1/gmm.V_1*gmm.B );
  m_xa = V_xa * ( gmm.A'*1/gmm.V_0*(y_i - gmm.a) + v_inv_divi * m_divi );
  m_xb = V_xb * ( gmm.B'*1/gmm.V_1*(y_i - gmm.b) + v_inv_divi * m_divi );
  V_ya = gmm.V_0 + gmm.A*v_divi*gmm.A';
  V_yb = gmm.V_1 + gmm.B*v_divi*gmm.B';
  m_ya = gmm.A*m_divi + gmm.a;
  m_yb = gmm.B*m_divi + gmm.b;

  % check for unnormalizable weights
  if (V_ya < 0) || (V_yb < 0)    
    err = MException('UpdateAugmentedGMM:MixWts', ...
        'Resulting mixture weights are not finite due to negative variance.');
    throw(err);
  end

  % compute association probability
  log_p = zeros(1,2);
  log_p(1) = log(gmm.w) - 0.5 * sign(V_ya) * log( 2*pi*abs(V_ya) ) - ...
    0.5 * ( y_i - m_ya )^2 * (1/V_ya);
  log_p(2) = log(1-gmm.w) - 0.5 * sign(V_yb) * log( 2*pi*abs(V_yb) ) - ...
    0.5 * ( y_i - m_yb )^2 * (1/V_yb);
  max_val = max( log_p );
  log_p = log_p - max_val;
  p = exp( log_p ) ./ sum( exp( log_p ) );
  
  % project moments
  [ V_aug, m_aug ] = ...
    gmm_moments( V_xa, V_xb, m_xa, m_xb, p );  
  
end

