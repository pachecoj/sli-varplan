function [ V_inv, m ] = gmm_moments( ...
  V_xa, V_xb, ...   % Component variances
  m_xa, m_xb, ...   % Component means
  p ...             % Normalized mixture weights
  )
% GMM_MOMENTS - Computes the mean and inverse variance of a two-component
%   Gaussian mixture.
%
% Jason L. Pacheco
% 12/15/11
%

  m = p(1) * m_xa + p(2) * m_xb;
  V = p(1) * V_xa + p(2) * V_xb + p(1) * (m - m_xa)*(m - m_xa)' + ...
    p(2) * (m - m_xb) * (m - m_xb)';
  V_inv = inv(V);    
  
end

