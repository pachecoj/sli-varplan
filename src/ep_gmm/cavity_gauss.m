function [ v_inv_divi, m_divi ] = cavity_gauss( ...
  v_inv_x, ...  % Posterior inverse variance
  m_x, ...      % Posterior mean
  lam_i, ...    % Message inverse variance
  eta_i ...     % Message scale
  )
% CAVITY_GAUSS - Compute cavity distribution for Gaussian EP.  Given
%   posterior moment parameters, and factor approximation canonical
%   parameters, returns moment parameters of the cavity distribution.
%
% Jason L. Pacheco
% 12/15/11
%

  v_inv_divi = v_inv_x - lam_i ;
  v_divi = 1/v_inv_divi;
  m_divi = v_divi * ( v_inv_x * m_x - eta_i );

end

