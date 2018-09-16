function pdf = dense_mixture_2D( ...
  gmm, ...
  x_vals, ...
  y_vals, ...  
  k ...
  )
  % DENSE_MIXTURE_2D - Create mixture PDF on dense grid

  %
  % FIRST COMPONENT
  %
  dy = y_vals(2) - y_vals(1);
  m = gmm.A*x_vals + gmm.a;
  D = pdist2( y_vals', m' );
  log_pdf_0 = -0.5 * log(2*pi*gmm.V_0) - 0.5 * D.^2 / gmm.V_0;
  log_pdf_0 = bsxfun(@minus, log_pdf_0, max(log_pdf_0,[],1));
  pdf_0 = bsxfun(@rdivide, exp(log_pdf_0), sum(exp(log_pdf_0),1)) ./ dy;
  
  %
  % SECOND COMPONENT
  %
  dist = abs(x_vals - gmm.loc(k));
  std_tgt = sqrt( gmm.c * dist + gmm.V_1 );  
  m = gmm.B*x_vals + gmm.b;
  D = pdist2( y_vals', m' );
  D = bsxfun(@rdivide, D, std_tgt);
  log_pdf_1 = bsxfun(@plus, - 0.5 * D.^2, -0.5*log(2*pi*(std_tgt.*std_tgt)) );
  log_pdf_1 = bsxfun(@minus, log_pdf_1, max(log_pdf_1,[],1));
  pdf_1 = bsxfun(@rdivide, exp(log_pdf_1), sum(exp(log_pdf_1),1)) ./ dy;
  
  %
  % MIXTURE
  %
  pdf = gmm.w * pdf_0 + (1-gmm.w) * pdf_1;
  
end