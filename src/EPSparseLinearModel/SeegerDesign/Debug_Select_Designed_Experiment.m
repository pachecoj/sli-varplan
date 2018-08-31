function [scores, scoresVI] = Debug_Select_Designed_Experiment(Nruns,num_Asamples,i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)

  n = size(X,1);
  
  % find all that candU that are not used yet
  poss = true(size(candU,2),1);
%   poss(order(1:i-1)) = false;
  poss = find(poss);
  poss_candU = candU(:,poss);  
  
  X = X';
  U = U';
  scores = zeros(size(poss_candU,2),Nruns);
  scoresVI = zeros(size(poss_candU,2),Nruns);    
  for irun = 1:Nruns
    xsamp = zeros(n, size(poss_candU,2), num_Asamples);
    asamp = zeros(n, n, num_Asamples);  
    
    for i = 1:num_Asamples
      % sample x_cand matix A from posterior
      [temp_candX, temp_candA] = eplin_sampxcand(X,U,sigma_noise^2,1,sitepi,siteb,L,gamma,poss_candU);
      xsamp(:,:,i) = temp_candX;
      asamp(:,:,i) = temp_candA;
      
      % compute score for
      for j = 1:n
        scores(:,irun) = scores(:,irun) + eplin_compscores(X,U(:,j),sigma_noise^2,1,...
          sitepi{j},siteb{j},L{j},gamma{j},temp_candX,poss_candU(j,:)');
      end          
    end    
    
    % variational bound w/ MCMC variance
    for ui=1:size(poss_candU,2)
      for j = 1:n
        
        % compute conditional precision
        xasamp = [ squeeze( xsamp(:, ui, :) ); squeeze(asamp(j,:,:)) ];
        m_xa = mean( xasamp, 2 );
        xx = bsxfun(@times, reshape(xasamp,2*n,1,[]), reshape(xasamp,1,2*n,[]));
        Exx = mean( xx, 3 );
        V_xa = Exx - m_xa * m_xa';        
        V_x = V_xa(1:n,1:n);
%         V_a = V_xa((n+1):end,(n+1):end);
%         Lam_xa = V_xa \ eye(2*n);
%         Lam_a_cond_x = Lam_xa((n+1):end,(n+1):end);
        
        % compute MI bound        
        Lam_a = L{j} * L{j}' / sigma_noise^2;
        scoresVI(ui,irun) = scoresVI(ui,irun) ...
          - 0.5 * sum( log( eig( Lam_a ) ) ) ...
          + 0.5 * sum( log( eig( V_x ) ) ) ...
          - 0.5 * sum( log( eig( V_xa ) ) );
      end
    end
  end
  scores = scores / n / num_Asamples;