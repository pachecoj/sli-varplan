function H = LLhessXA(x, A, u, m, Lam, sigma_noise)
  n = length(m);
  Hxx = zeros(n,n);
  Haa = zeros(n^2,n^2);
  Hxa = zeros(n,n^2);
    
  for i=1:n
    a = A(i,:)';
    idx_start = 1 + n*(i-1);
    idx_end = idx_start + n - 1;
    Hxx = Hxx - a * a' / sigma_noise^2;
    Haa(idx_start:idx_end,idx_start:idx_end) = ...
      - Lam{i} - x*x'/sigma_noise^2;
    
    % cross terms
    H_xa_blk = - a*x'/sigma_noise^2;
    H_xa_diag = u(i)/sigma_noise^2 - a'*x/sigma_noise^2;
    Hxa(:,idx_start:idx_end) = H_xa_blk + H_xa_diag * eye(n);
         
  end
  
  % build hessian blockwise
  H = blkdiag(Hxx, Haa);
  H(1:n,(n+1):end) = Hxa;
  H((n+1):end,1:n) = Hxa';
end