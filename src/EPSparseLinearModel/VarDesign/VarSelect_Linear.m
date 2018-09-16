function [ind,scores] = VarSelect_Linear(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
if i <= num_initial_rand
  ind = Select_Random_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order);
  scores = 0;
else
  n = size(X,1);

  % find all the candU that are not used yet
  poss = true(size(candU,2),1);
  poss(order(1:i-1)) = false;
  poss = find(poss);
  poss_candU = candU(:,poss);
  X = X';
  U = U';
  
  % sample matrix A from posterior
  N_samp = 100;
  A = zeros([n,n,N_samp]);  
  for k = 1:N_samp    
    [~, thisA] = eplin_sampxcand(X,U,sigma_noise^2,1,sitepi,siteb,L,gamma,poss_candU);    
    A(:,:,k) = thisA;
  end
  
  % compute score for each U candidate
  scores = zeros(size(poss_candU,2),1);  
  for ui = 1:size(poss_candU,2)
    u = poss_candU(:,ui);
    
    % compute moments
    Gamma = zeros(n,n);    
    M = zeros(n,N_samp);
    for k=1:N_samp    
      Gamma = Gamma + A(:,:,k)\(sigma_noise^2 + u*u')/(A(:,:,k)');    
      M(:,k) = A(:,:,k)\u;
    end
    mbar = mean( M, 2 );
    S = 1/N_samp * Gamma - mbar * mbar';
            
    % sum score over rows of A
    for j = 1:n
      calA = shiftdim( A(j,:,:), 1 ); % NxM      
      AMT = calA * M';
            
      % compute parameters
      Lam_n = L{j} * L{j}' / sigma_noise / sigma_noise;
      h_n = L{j}' \ gamma{j};
      F_n = Lam_n * ( 1/N_samp * AMT - h_n * mbar' ) / S;
      b_n = Lam_n * h_n - F_n * mbar;

      % compute score
      scores(ui) = scores(ui) + 1/N_samp * trace(F_n * AMT') + h_n' * b_n;
      scores(ui) = scores(ui) - 1/2 * 1/N_samp * trace( F_n' / Lam_n * F_n * Gamma );
      scores(ui) = scores(ui) - mbar' * F_n' / Lam_n * b_n;
      scores(ui) = scores(ui) - 1/2 * b_n' / Lam_n * b_n;    
      
    end
  end
%   scores = scores / n / num_Asamples;
  [void,ind] = max(scores);
  ind = poss(ind);
  fprintf('Selected candidate %i with score %g.\n',ind,void);
end


