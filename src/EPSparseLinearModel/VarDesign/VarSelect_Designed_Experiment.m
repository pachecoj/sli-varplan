function ind = VarSelect_Designed_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
if i <= num_initial_rand
  ind = Select_Random_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order);
else
  n = size(X,1);

  % find all that candU that are not used yet
  poss = true(size(candU,2),1);
  poss(order(1:i-1)) = false;
  poss = find(poss);
  poss_candU = candU(:,poss);
  X = X';
  U = U';
  
  % sample matrix A from posterior
  num_Asamples = 40;
  A = zeros([n,n,num_Asamples]);  
  Ainv = zeros([n,n,num_Asamples]);  
  AinvT = zeros([n,n,num_Asamples]);    
  for k = 1:num_Asamples    
    [~, thisA] = eplin_sampxcand(X,U,sigma_noise^2,1,sitepi,siteb,L,gamma,poss_candU);    
    A(:,:,k) = thisA;
    Ainv(:,:,k) = inv( thisA );
    AinvT(:,:,k) = Ainv(:,:,k)';
  end
  
  % compute score for each U candidate
  scores = zeros(size(poss_candU,2),1);
  for ui = 1:size(poss_candU,2)
    u = poss_candU(:,ui);
    
    % compute moments
    Lam = zeros(n,n);    
    M = zeros(n,num_Asamples);
    for k=1:num_Asamples
      Lam = Lam + sigma_noise^2*Ainv(:,:,k)*AinvT(:,:,k) + ...
        Ainv(:,:,k)*(u*u')*AinvT(:,:,k);
      M(:,k) = Ainv(:,:,k)*u;      
    end
    invLam = inv(Lam);
            
    % sum score over rows of A
    for j = 1:n
      calA = shiftdim( A(j,:,:), 1 ); % NxM      
      AM = calA * M';
      scores(ui) = scores(ui) + trace( AM * invLam * AM' );
    end
  end
%   scores = scores / n / num_Asamples;
  [void,ind] = max(scores);
  ind = poss(ind);
  fprintf('Selected candidate %i with score %g.\n',ind,void);
end


