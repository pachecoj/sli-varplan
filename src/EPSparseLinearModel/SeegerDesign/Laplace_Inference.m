function [Qm,Qv,times,order] = Laplace_Inference(U,X,candU,candX,sigma_noise,tau,num_inclus,orderfunction)
n = size(X,1);

% initialise n sparse linear models
sitepi = cell(n,1);
siteb = cell(n,1);
L = cell(n,1);
gamma = cell(n,1);
for i = 1:n
  sitepi{i} = repmat(0.5*(tau*sigma_noise)^2,n,1);
  siteb{i} = zeros(n,1);
  L{i} = zeros(n);
  gamma{i} = zeros(n,1);
end
pithres = 1e-16;
wkvec1 = zeros(n,1);wkvec2 = zeros(n,1);wkvec3 = zeros(n,1);
Qm = cell(num_inclus+1,1);
Qv = cell(num_inclus+1,1);

% do the initial updates
fprintf('Initialisation ...\n');
% fprintf('WARNING: Doing ADF inference.\n');
for i = 1:n
  for j = 1:20
    sweep_ord = randperm(n);
    delta=ep_splinsweep(X',U(i,:),sigma_noise^2,tau*sigma_noise,...
      sitepi{i},siteb{i},sweep_ord,L{i},gamma{i},pithres,...
      double(j==1),wkvec1,wkvec2,wkvec3,0,1,0.5);
    %         For checking whether the algotithm converges
%             fprintf(1,'i = %i j = %i delta=%f\n',i,j,delta);
  end
end
[Qm{1},Qv{1}] = eval_marginal_mean_var(L,gamma,sigma_noise);

% now iterate over experiments
order = zeros(num_inclus,1);
for m = 1:num_inclus
  fprintf('Performing %i-th measurement (Design + Update) ...\n',m);
  [order(m), scores] = orderfunction(m,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order);
  
%   % DEBUG
%   if m>20
%     tmp_order = order;
%     tmp_order(m) = 0;
%     [~, scores2] = Select_Designed_Experiment(m,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,tmp_order,20);
%     [~,idx] = sort( scores2, 'descend' );
%     figure;
%     hold on;
%     plot(scores(idx))    
%     keyboard;
%   end  
  
  X = [X,candX(:,order(m))];
  U = [U,candU(:,order(m))];
%   fprintf('WARNING: Doing ADF inference.\n');
  for i = 1:n
    for j = 1:20
      sweep_ord = randperm(n);
      delta=ep_splinsweep(X',U(i,:),sigma_noise^2,tau*sigma_noise,...
        sitepi{i},siteb{i},sweep_ord,L{i},gamma{i},pithres,...
        double(j==1),wkvec1,wkvec2,wkvec3,0,1,0.5);
      %         For checking whether the algotithm converges
%                   fprintf(1,'i = %i j = %i delta=%f\n',i,j,delta);
    end
  end
  [Qm{m+1},Qv{m+1}] = eval_marginal_mean_var(L,gamma,sigma_noise);
end
times = 0;