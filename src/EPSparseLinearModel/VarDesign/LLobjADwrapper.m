function [negLogP, g, H] = LLobjADwrapper(x0, u, m, Lam, V, sigma_noise)
  n = numel( x0 );

  % evaluate objective
  negLogP = LLobjAD(x0, u, m, Lam(:), V(:), sigma_noise);
  
  % evaluate gradient
  g = []; H = [];
  if nargout > 1
    x.f = x0;
    x.dx = ones( size( x0 ) );
    gadi = LLobjAD_ADiGatorGrd(x, u, m, Lam(:), V(:), sigma_noise);
    g = gadi.dx;
    
    if nargout > 2
      Hadi = LLobjAD_ADiGatorHes(x, u, m, Lam(:), V(:), sigma_noise);
      H = reshape( Hadi.dxdx, [n, n] );
    end
  end  
end