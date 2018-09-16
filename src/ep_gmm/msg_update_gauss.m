function [lam_i_new, eta_i_new] = msg_update_gauss(...
    v_inv_x, ...    % posterior inverse variance
    m_x, ...        % posterior mean
    v_inv_divi, ... % cavity inverse variance
    m_divi ...      % cavity mean
  )
% MSG_UPDATE_GAUSS - Updates Gaussian message given posterior moment
%   parameters and canonical moment parameters.  Returns message canonical
%   parameters.
%
% RETURNS:
%   lam_i_new - Message inverse variance (precision) parameter
%   eta_i_new - Message length scale parameter
%
% Jason L. Pacheco
% 12/14/11
%
 
  if ( abs( v_inv_x - v_inv_divi ) > eps )
    lam_i_new = (v_inv_x - v_inv_divi);
    eta_i_new = v_inv_x * m_x - v_inv_divi * m_divi;
  else
    lam_i_new = 0;
    eta_i_new = 0;
  end
  
end