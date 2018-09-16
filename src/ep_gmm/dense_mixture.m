function pdf = dense_mixture( ...
  gmm, ...
  x_vals, ...
  y, ...  
  k ...
  )
  % DENSE_MIXTURE - Create mixture PDF on dense grid

  dist = abs(x_vals - gmm.loc(k));
  var_tgt = gmm.c * dist + gmm.V_1;
  pdf = gmm.w .* normpdf(y, gmm.A*x_vals + gmm.a, sqrt(gmm.V_0))' + ...
    (1-gmm.w) .* normpdf(y, gmm.B*x_vals + gmm.b, sqrt(var_tgt))' ;
end