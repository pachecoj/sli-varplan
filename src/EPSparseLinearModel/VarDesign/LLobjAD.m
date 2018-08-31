function negLogP = LLobjAD(x, u, m, Lam, V, sigma_noise)

  n = length(m);
  Lam = reshape( Lam, [ n, n * n ] );
  V = reshape( V, [ n, n * n ] );
  
  % EP posterior
  logP = 0;
  A = zeros(n,n);
  colstart = 1;
  for i = 1:n
%     t_start = tic();    
    thisLam = Lam(:, colstart:(colstart + n - 1));
    thisV = V(:, colstart:(colstart + n - 1));
    this_m = m(:,i);    
    this_h = thisLam * this_m;        
    a = ( thisLam + (x * x') / sigma_noise^2 ) \ ( this_h + x*u(i)/sigma_noise^2);
%     invLamXX = thisV - (1/sigma_noise^2) * thisV * (x * x') * thisV / (1 + x' * thisV * x / sigma_noise^2 );
%     a = invLamXX * ( this_h + x*u(i)/sigma_noise^2);
    logP = logP - 0.5 * (a - this_m)' * thisLam * (a - this_m); % + -n/2 * log(2*pi) + 0.5 * log(det(Lam{i})) ...      
    A(i,:) = a;
%     t_stop = toc(t_start);
%     fprintf('%d (%0.3fs)\n', i, t_stop);
    colstart = colstart + n;
  end
  
  % data term
  logP = logP - 0.5/sigma_noise/sigma_noise * (u-A*x)' * (u-A*x); % -n/2 * log(2*pi*sigma_noise*sigma_noise) ...    
  negLogP = - logP;
end