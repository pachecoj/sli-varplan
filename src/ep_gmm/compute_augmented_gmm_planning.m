function [ V_a, V_b, m_a, m_b ] = compute_augmented_gmm_planning( ...
  gmm, ...   % Gauss. Mixture Model params
  invV, ...  % Posterior inverse variance of X
  eta ...    % Posterior canonical location-scale of X  
  )
% COMPUTE_AUGMENTED_GMM_PLANNING - Multiplies Gaussian posterior p(x) with
%   two-component GMM likelihood to form joint pdf p(X,Y).
%
% J. Pacheco, 2017
%

  % unpack stuff
  [ A, B, a, b, w, V_0, V_1 ] = deal( gmm.A, gmm.B, gmm.a, gmm.b, ...
    gmm.w, gmm.V_0, gmm.V_1 );
  invV_0 = 1/V_0;
  invV_1 = 1/V_1;

  %
  % FIRST COMPONENT
  %
  invV_a = [ ...
    A^2 * invV_0 + invV,  A*invV_0; ...
    A*invV_0, invV_0 ...
    ];
  V_a = inv( invV_a );
  m_a = V_a * [ eta; 0 ];
  
  %
  % SECOND COMPONENT
  %
  invV_b = [ ...
    B^2 * invV_1 + invV,  B*invV_1; ...
    B*invV_1,  invV_1 ...
    ];
  V_b = inv( invV_b );
  m_b = V_b * [ eta; 0 ];

end