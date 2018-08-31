function ind = Select_Random_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order)
  m = size(candU,2);n = size(candU,1);
  % find all that candU that are not used yet
  poss = true(m,1);
%   poss(order(1:i-1)) = false;
  poss = find(poss);
  % choose a random one
  ind = irand(1,1,length(poss));
  ind = poss(ind);
  
function res = irand(n,min,max)
  res = floor((rand(n,1)-1e-10)*(max-min)) + min;   