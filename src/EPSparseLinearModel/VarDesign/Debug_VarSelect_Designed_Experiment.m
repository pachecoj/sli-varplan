function [scores, scoresLD] = Debug_VarSelect_Designed_Experiment(Nruns, i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
  n = size(X,1);

  % find all that candU that are not used yet
  poss = true(size(candU,2),1);
  poss(order(1:i-1)) = false;
  poss = find(poss);
  poss_candU = candU(:,poss);
  X = X';
  U = U';
  
  % sample matrix A from posterior
  num_Asamples = 500;
  scores = zeros(size(poss_candU,2),Nruns);  
  scoresLD = zeros(size(poss_candU,2),Nruns);  
  for irun = 1:Nruns
    fprintf('Run: %d\n',irun);
    A = zeros([n,n,num_Asamples]);
    Ainv = zeros([n,n,num_Asamples]);
    AinvT = zeros([n,n,num_Asamples]);
    for k = 1:num_Asamples
      [temp_candX, thisA] = eplin_sampxcand(X,U,sigma_noise^2,1,sitepi,siteb,L,gamma,poss_candU);
      A(:,:,k) = thisA;
      Ainv(:,:,k) = inv( thisA );
      AinvT(:,:,k) = Ainv(:,:,k)';
      
      % compute score using Seeger design
      for j = 1:n
        scoresLD(:,irun) = scoresLD(:,irun) + eplin_compscores(X,U(:,j),sigma_noise^2,1,...
          sitepi{j},siteb{j},L{j},gamma{j},temp_candX,poss_candU(j,:)');
      end
    end
    
    % compute score for each U candidate    
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
        scores(ui,irun) = scores(ui,irun) + trace( AM * invLam * AM' );
      end
    end
    
    
  end
%   scores = scores / n / num_Asamples;