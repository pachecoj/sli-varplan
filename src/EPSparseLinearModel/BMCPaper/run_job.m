function run_job(dir,run_opt,rounds)

% some overhead, making sure that cache files are deleted
if ~exist(dir,'dir'), mkdir(dir), end;
evalfile = sprintf('%s/%s-evaluation.mat',dir,get_run_opt_FN(run_opt));
if exist(evalfile), delete(evalfile), end;

N = run_opt.data_opt.N;
for r = 1:rounds
  fprintf('\n\nRunning round %i ... \n',r);
  load(sprintf('%s/%s-r=%i-data.mat',dir,run_opt.data_opt.name,r),'A','U','X','candU','candX','XN');
  
  % estimate sigma from XN or use it as given
%   if 0 < run_opt.sigma_method && run_opt.sigma_method < 1
    sigma_noise = run_opt.sigma_method;
    fprintf('Sigma directly given: sigma = %g\n',sigma_noise);
%   else
%     switch run_opt.sigma_method
%       case -1 % known (not necessarily true for SDE simluation)
%         sigma_noise = run_opt.data_opt.sigma_noise;
%       case 1
%         sigma_noise = std(XN(:));
%       case 2
%         sigma_noise = mean(std(XN,0,2) ./ abs(diag(A)));
%       case 3 % old and working
%         sigma_noise = std(XN(:)) ./ mean(abs(diag(A)));
%       case 4
%         x = log10(run_opt.data_opt.input_scale);
%         y = -9 +  7/4* (x+4);
%         sigma_noise = 10^(y);
%       case 5
%         ee = A * [X,candX] - [U,candU];
%         sigma_noise = std(ee(:));
%       otherwise
%         error('Unkown method to determine Sigma.');
%     end
%     fprintf('Sigma method %i: sigma = %g\n',run_opt.sigma_method,sigma_noise);
%   end
  
  % estimate tau or sigma respectively
  GG = (A > run_opt.eps) | (A < -run_opt.eps);
  deg = mean(sum(GG,2));
  switch run_opt.method
    case {1,3,6,7}
      if run_opt.tau ~= -1
        tau_est = run_opt.tau;
      else
        tau_est = -1/run_opt.eps * log(deg/N);
        fprintf('Estimated Laplace sigma: %g (tau = %g)\n',sqrt(2)/tau_est,tau_est);
      end
    case {4,5}
      if run_opt.tau ~= -1
        sigma_model = 1/run_opt.tau;
      else
        opt = optimset('largescale','off','Display','off');
        sigma_model = fminunc(@(x)((1-mycdf(0,x,run_opt.eps)+mycdf(0,x,-run_opt.eps)-deg/N)^2),1,opt);
        fprintf('Estimated Gaussian sigma: %g\n',sigma_model);
      end
  end
  
  tic
  switch run_opt.method
    case 1
      [Qm,Qv,times,order] = Laplace_Design(U,X,candU,candX,sigma_noise,tau_est,run_opt.num_inclus,run_opt.initial_rand);
    case 3
      [Qm,Qv,times,order] = Laplace_Random(U,X,candU,candX,sigma_noise,tau_est,run_opt.num_inclus);
    case 4
      [Qm,Qv,times,order] = Gauss_Design(U,X,candU,candX,sigma_noise,sigma_model,run_opt.num_inclus,run_opt.initial_rand);
    case 5
      [Qm,Qv,times,order] = Gauss_Random(U,X,candU,candX,sigma_noise,sigma_model,run_opt.num_inclus);
    case 6
      [Qm,Qv,times,order] = Laplace_VarDesign_Affine(U,X,candU,candX,sigma_noise,tau_est,run_opt.num_inclus,run_opt.initial_rand);
    case 7
      [Qm,Qv,times,order] = Laplace_VarDesign_Linear(U,X,candU,candX,sigma_noise,tau_est,run_opt.num_inclus,run_opt.initial_rand);
    otherwise
      error('Unkown method');
  end
  times = toc;
  fprintf('Used %g sec for computing experiment.\n',times);
  
  round_fn = sprintf('%s/%s-r=%i-results.mat',dir,get_run_opt_FN(run_opt),r);
  save(round_fn, 'Qm','Qv','times','order');
  clear Qm Qv times order U X candU candX
end


function [Qm,Qv,times,order] = Laplace_Design(U,X,candU,candX,sigma_noise,tau,num_inclus,num_initial_rand)
orderfunction = @(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order)...
  (Select_Designed_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand));
[Qm,Qv,times,order] = Laplace_Inference(U,X,candU,candX,sigma_noise,tau,num_inclus,orderfunction);

function [Qm,Qv,times,order] = Laplace_VarDesign_Affine(U,X,candU,candX,sigma_noise,tau,num_inclus,num_initial_rand)
orderfunction = @(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order)...
  (VarSelect_Affine(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand));
[Qm,Qv,times,order] = Laplace_Inference(U,X,candU,candX,sigma_noise,tau,num_inclus,orderfunction);

function [Qm,Qv,times,order] = Laplace_VarDesign_Linear(U,X,candU,candX,sigma_noise,tau,num_inclus,num_initial_rand)
orderfunction = @(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order)...
  (VarSelect_Linear(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand));
[Qm,Qv,times,order] = Laplace_Inference(U,X,candU,candX,sigma_noise,tau,num_inclus,orderfunction);

function [Qm,Qv,times,order] = Laplace_Random(U,X,candU,candX,sigma_noise,tau,num_inclus)
[Qm,Qv,times,order] = Laplace_Inference(U,X,candU,candX,sigma_noise,tau,num_inclus,@Select_Random_Experiment);

function [Qm,Qv,times,order] = Gauss_Design(U,X,candU,candX,sigma_noise,sigma_model,num_inclus,num_initial_rand)
orderfunction = @(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order)...
  (Select_Designed_Experiment(i,X,U,sitepi,siteb,L,gamma,candU,sigma_noise,order,num_initial_rand));
[Qm,Qv,times,order] = Gauss_Inference(U,X,candU,candX,sigma_noise,sigma_model,num_inclus,orderfunction);

function [Qm,Qv,times,order] = Gauss_Random(U,X,candU,candX,sigma_noise,sigma_model,num_inclus)
[Qm,Qv,times,order] = Gauss_Inference(U,X,candU,candX,sigma_noise,sigma_model,num_inclus,@Select_Random_Experiment);










