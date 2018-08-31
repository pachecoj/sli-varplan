function g = LLgradXA(x, A, u, m, Lam, sigma_noise)
  n = length(m);    
  gx = 0;
  gA = zeros(size(A));
  
  for i = 1:n    
    a = A(i,:)';
    gx = gx + a*u(i) - a*x'*a ;
    gA(i,:) = -( Lam{i} + x*x'/sigma_noise^2 )*a + x*u(i)/sigma_noise^2 + Lam{i}*m{i};
%     gA(i,:) = - Lam{i}*a - x*a'*x/sigma_noise^2 + x*u(i)/sigma_noise^2 + Lam{i}*m{i};
  end
  gx = gx/sigma_noise/sigma_noise;
  
  % repack gradient
  gAt = gA';
  g = [ gx(:); gAt(:) ];
end