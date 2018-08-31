function A = psdproj( A, min_eigval )
  A = 0.5 * (A + A'); % make symmetric
  [Q,D] = eig( A );  
  evals = diag( D );
  min_eigval = 1e-15 * max(evals);
  evals( evals < min_eigval ) = min_eigval;
  Dproj = diag( evals );
  A = Q*Dproj*Q';
  A = 0.5 * (A + A'); % make symmetric
end