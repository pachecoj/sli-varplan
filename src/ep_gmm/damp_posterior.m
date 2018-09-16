function [ v_inv_damp, m_damp ] = damp_posterior( ...
  v_inv_new, ...  % (new) inverse variance
  m_new, ...      % (new) mean
  v_inv_old, ...  % (old) inverse variance
  m_old, ...      % (old) mean
  alpha ...       % stepsize
  )
% DAMP_POSTERIOR Applies linear update damping to Gaussian canonical
%   parameters corresponding to the specified mean parameters.  That is,
%   this function converts the input mean parameters to canonical
%   parameters of a Gaussian, then returns a linear combination of the
%   (new) and (old) canonical parameters.
% 
% Jason L. Pacheco
% 12/15/11
%

  % check if stepsize=1
  if alpha == 1.0
    v_inv_damp = v_inv_new;
    m_damp = m_new;
    return;
  end

  % otherwise, damp
  eta_damp = alpha * v_inv_new*m_new + (1-alpha)*v_inv_old*m_old;
  lambda_damp = alpha * v_inv_new + (1-alpha)*v_inv_old;
  
  % convert to mean params
  v_inv_damp = lambda_damp;
  m_damp = (1/lambda_damp) * eta_damp;

end

