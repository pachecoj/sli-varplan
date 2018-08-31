function [Qm,Qv,times,order] = Gauss_Inference(U,X,candU,candX,sigma_noise,sigma_model,num_inclus,orderfunction)
n = size(X,1);

% initialise n sparse linear models
sitepi = cell(n,1);
siteb = cell(n,1);
L = cell(n,1);
gamma = cell(n,1);
for i = 1:n
  sitepi{i} = repmat((sigma_noise/sigma_model)^2,n,1);
  siteb{i} = zeros(n,1);
  L{i} = zeros(n);
  gamma{i} = zeros(n,1);
end
Qm = cell(num_inclus+1,1);
Qv = cell(num_inclus+1,1);

% do the initial updates
fprintf('Initialisation ...\n');
for i = 1:n
  L{i} = chol(X*X' + diag(sitepi{i}))';
  gamma{i} = L{i} \ (X*U(i,:)');
end
[Qm{1},Qv{1}] = eval_marginal_mean_var(L,gamma,sigma_noise);

% now iterate over experiments
order = zeros(num_inclus,1);
for m = 1:num_inclus
  fprintf('Performing %i-th measurement (Design + Update) ...\n',m);
  order(m) = orderfunction(m,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order);
  
  X = [X,candX(:,order(m))];
  U = [U,candU(:,order(m))];
  for i = 1:n
    L{i} = chol(X*X' + diag(sitepi{i}))';
    gamma{i} = L{i} \ (X*U(i,:)');
  end
  [Qm{m+1},Qv{m+1}] = eval_marginal_mean_var(L,gamma,sigma_noise);
end
times = 0;