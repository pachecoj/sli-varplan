function ind = Select_Designed_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand)
if i <= num_initial_rand
  ind = Select_Random_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order);
else
  n = size(X,1);
  % find all that candU that are not used yet
  poss = true(size(candU,2),1);
%   poss(order(1:i-1)) = false;
  poss = find(poss);
  poss_candU = candU(:,poss);
  
  
  X = X';
  U = U';
  scores = zeros(size(poss_candU,2),1);
  num_Asamples = 20; 
  for i = 1:num_Asamples
    % sample x_cand matix A from posterior
    temp_candX = eplin_sampxcand(X,U,sigma_noise^2,1,sitepi,siteb,L,gamma,poss_candU);
    
    % compute score for
    for j = 1:n
      scores = scores + eplin_compscores(X,U(:,j),sigma_noise^2,1,...
        sitepi{j},siteb{j},L{j},gamma{j},temp_candX,poss_candU(j,:)');
    end
  end
  scores = scores / n / num_Asamples;
  [void,ind] = max(scores);
  ind = poss(ind);
  fprintf('Selected candidate %i with score %g.\n',ind,void);
end