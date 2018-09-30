% DEMO_VARPLAN - Demonstrate variational inference & planning for simple
%   Gaussian mixture model.
%
% At each time scan data are sampled from a two-component Gaussian mixture
% model with possible K sensors.  A latent parameter is first sampled:
%
%  x ~ N( m0_x, V0_x)
%
% At each time closed-loop greedy planning selects one of K sensors and
% draws an observation from a two-component Gaussian mixture:
%
%  y | x; k ~ w * N( Ax + a, V_0) + (1-w) * N( Bx + b, V_1k)
%
% Where V_0 is a fixed noise variance and V_1k is a target observation
% variance which depends on the distance between target and sensor k:
%
% V_1k = c * abs( loc(k) - x ) + V_1
%
% Here loc(k) is the location on sensor k (fixed) and c is a scaling
% parameter which is the same for all sensors.  This script simulates 
% multiple trials and compares inference and planning under various 
% approximations ( ADF / EP / MCMC ).
%
% J. Pacheco, 2017
%

clear variables
close all

seed = 1;
rng(seed);

% INIT. PARAMETERS
D = 1;
A = 0;
B = 1;
a = 0;
b = 0;
w = 0.25;
V_0 = 10;
V_1 = 1;
c = 1;
m0_x = 0;
V0_x = 100;
num_data = 20;
Nsensors = 10;
sensorBounds = [-25, 25];
model = init_gmm(A, B, a, b, w, c, V_0, V_1, m0_x, V0_x, Nsensors, sensorBounds);
N = num_data;

% INIT. CONTROL FLAGS
plot_flag = true;
max_iters = 100;
conv_thresh = 1E-6;
damp_type = 1;
min_damp_fact = 0.25;
ep_opt = init_ep_opt(...
  1:N, max_iters, conv_thresh, damp_type, min_damp_fact, plot_flag);
adf_opt = ep_opt;
adf_opt.max_iters = 1;
Nsamp = 50;

% INIT. OTHER STUFF
x_grid = linspace(-50, 50, 1000 );
dx = x_grid(2) - x_grid(1);

%
% DO MULTIPLE TRIALS
%
Nruns = 1;
MI_err_mcmc = zeros(Nsensors,N,Nruns);
H_cond_err_ep = zeros(Nsensors,N,Nruns);  
H_cond_err_adf = zeros(Nsensors,N,Nruns);
err_mcmc = zeros(N,Nruns);
err_ep = zeros(N,Nruns);
err_adf = zeros(N,Nruns);
err_mcmc_state = zeros(N,Nruns);
err_ep_state = zeros(N,Nruns);
err_adf_state = zeros(N,Nruns);
err_random_state = zeros(N,Nruns);
for irun = 1:Nruns
  fprintf('Run %d of %d.\n', irun, Nruns);
    
  % SAMPLE DATA  
  [x, y_all] = sample_gmm(model, N);  
  
  % CREATE FIGURE HANDLES
  if ((ep_opt.plot_flag) && (Nruns == 1))
    h_fig = figure('InvertHardcopy','off','Color',[1 1 1]);
    set(h_fig,'Units','Pixels');
    set(h_fig, 'Position', [0 0 800 600]);
    %   h_planfig = figure('InvertHardcopy','off','Color',[1 1 1]);
  else
    h_fig = [];
  end    
  
  %
  % EXPECTATION PROPAGATION
  %
  k_vec_ep = 1;
  for n=1:(N-1)
    idx = sub2ind(size(y_all),k_vec_ep,1:n);
    this_y = y_all(idx);
    ep_opt.sched = 1:n;
    
    % calculate true posterior
    [true_posterior_pdf, ~, best_fit_pdf] = true_posterior_gmm( model, this_y, x_grid, k_vec_ep );
    
    % inference
    fnplot = @(h_fig, y_hat_pdf, msg_pdf, y, x_vals) plot_ep_iter( ...
      h_fig, true_posterior_pdf, y_hat_pdf, best_fit_pdf, msg_pdf, [], y, x_vals );
    ep_state = ep_gmm(model, ep_opt, x, this_y, k_vec_ep, x_grid, fnplot, h_fig);
    err_ep_state(n,irun) = abs( ep_state.m_x - x );
    if (ep_opt.max_iters>1)
      if ep_state.iters < ep_opt.max_iters
        fprintf('EP Converged in %d iters.\n', ep_state.iters);
      else
        fprintf('EP did not converge.\n');
      end
    end
    
    % planning
    [ k_vec_ep(n+1), H_cond_err_ep(:,n+1,irun) ] = vi_plan_gmm( model, x_grid, ep_state, h_fig );
        
    % adjust plot
    if (ep_opt.plot_flag && (Nruns==1))
      subplot(2,1,1);
      set(gca,'FontSize',14);
      xlim([-25 25]);
      legend('True Posterior','ADF','Best Gaussian','Messages','Observations');
      title(sprintf('EP (%d Observations)',n));
      pause(1);
    end    
  end
  
  
  %
  % ASSUMED DENSITY FILTERING
  %  
  k_vec_adf = 1;  
  for n=1:(N-1)
    idx = sub2ind(size(y_all),k_vec_adf,1:n);
    this_y = y_all(idx);
    adf_opt.sched = 1:n;
    
    % calculate true posterior
    [true_posterior_pdf, ~, best_fit_pdf] = true_posterior_gmm( model, this_y, x_grid, k_vec_adf );
    
    % inference    
    fnplot = @(h_fig, y_hat_pdf, msg_pdf, y, x_vals) plot_ep_iter( ...
      h_fig, true_posterior_pdf, y_hat_pdf, best_fit_pdf, msg_pdf, [], y, x_vals );
    adf_state = ep_gmm(model, adf_opt, x, this_y, k_vec_adf, x_grid, fnplot, h_fig);
    err_adf_state(n,irun) = abs( adf_state.m_x - x );
    
    % planning
    [ k_vec_adf(n+1), H_cond_err_adf(:,n+1,irun) ] = vi_plan_gmm( model, x_grid, adf_state, h_fig );
    
    % adjust plot
    if (adf_opt.plot_flag && (Nruns==1))
      subplot(2,1,1);
      set(gca,'FontSize',14);
      xlim([-25 25]);
      legend('True Posterior','EP','Best Gaussian','Messages','Observations');
      title(sprintf('ADF (%d Observations)',n));
      pause(1);
    end    
  end
    
  
  %
  % MCMC
  %
  k_vec_mcmc = 1;  
  for n=1:(N-1)
    fprintf('MCMC %d of %d\n',n,N-1);
    idx = sub2ind(size(y_all),k_vec_mcmc,1:n);
    this_y = y_all(idx);
    
    % calculate true posterior    
    [true_posterior_pdf, ~, best_fit_pdf] = true_posterior_gmm(model, this_y, x_grid, k_vec_mcmc);
    
    % inference
    x_mcmc = mcmc_gmm(Nsamp, model, x, this_y, k_vec_mcmc, x_grid, [], h_fig);
    err_mcmc_state(n,irun) = abs( mean(x_mcmc) - x );
    
    % planning
    [ k_vec_mcmc(n+1), MI_err_mcmc(:,n+1,irun) ] = ...
      mcmc_plan_gmm(Nsamp, model, x_grid, x_mcmc, true_posterior_pdf, h_fig);
    
    % plot stuff    
    if (ep_opt.plot_flag && (Nruns==1))
      subplot(2,1,1);
      cla;
      set(gca,'FontSize',14);
      plot(x_grid, true_posterior_pdf, '-b', 'LineWidth', 2);
      hold on
      plot(x_mcmc, 0.01 * ones(numel(x_mcmc)), 'ok', 'MarkerSize', 10);
      xlim([-25 25]);
      legend('True Posterior','MCMC Samples');        
      title(sprintf('MCMC (%d Observations)',n));
      pause(1);
    end      
  end  
  
  %
  % RANDOM
  %
  k_vec_rand = 1;  
  for n=1:(N-1)
    fprintf('Random %d of %d\n',n,N-1);
    idx = sub2ind(size(y_all),k_vec_rand,1:n);
    this_y = y_all(idx);
    
    % calculate true posterior    
    [true_posterior_pdf, ~, best_fit_pdf] = true_posterior_gmm(model, this_y, x_grid, k_vec_rand);
    
    % inference
    x_rand = mcmc_gmm(Nsamp, model, x, this_y, k_vec_rand, x_grid, [], h_fig);
    err_random_state(n,irun) = abs( mean(x_rand) - x );
    
    % planning
    k_vec_rand(n+1) = randi(Nsensors);     
  end  
  
  % compute errors
  idx_mcmc = sub2ind( size( MI_err_mcmc ), k_vec_mcmc, 1:N, irun*ones(1,N) );
  err_mcmc(:,irun) = abs( MI_err_mcmc( idx_mcmc ) );
  idx_ep = sub2ind( size( H_cond_err_ep ), k_vec_ep, 1:N, irun*ones(1,N) );
  err_ep(:,irun) = abs( H_cond_err_ep( idx_ep ) );
  idx_adf = sub2ind( size( H_cond_err_adf ), k_vec_adf, 1:N, irun*ones(1,N) );
  err_adf(:,irun) = abs( H_cond_err_adf( idx_adf ) );
end

% plot MI errs (all)
figure('InvertHardcopy','off','Color',[1 1 1]);
set(gca,'FontSize',14);
hold on
for irun = 1:Nruns
  plot(2:N, err_mcmc(2:N,irun), '-b', 'LineWidth', 1);
  plot(2:N, err_ep(2:N,irun), '-r', 'LineWidth', 1);
  plot(2:N, err_adf(2:N,irun), '-k', 'LineWidth', 1);  
end
xlim([2 N])
legend('MCMC','EP','ADF');
xlabel('Observations');
ylabel('MI Estimate Error')

% plot MI errs (stats)
figure('InvertHardcopy','off','Color',[1 1 1]);
set(gca,'FontSize',14);
hold on
err_mcmc_mean = mean(err_mcmc(2:N,:), 2);
err_mcmc_std = std(err_mcmc(2:N,:), 0, 2);
err_ep_mean = mean(err_ep(2:N,:), 2);
err_ep_std = std(err_ep(2:N,:), 0, 2);
err_adf_mean = mean(err_adf(2:N,:), 2);
err_adf_std = std(err_adf(2:N,:), 0, 2);
plot(2:N, err_mcmc_mean, '-b', 'LineWidth', 2);
plot(2:N, err_ep_mean, '-r', 'LineWidth', 2);
plot(2:N, err_adf_mean, '-k', 'LineWidth', 2);  
plot(2:N, err_mcmc_mean - err_mcmc_std, '--b', 'LineWidth', 1);
plot(2:N, err_ep_mean - err_ep_std, '--r', 'LineWidth', 1);
plot(2:N, err_adf_mean - err_adf_std, '--k', 'LineWidth', 1);  
plot(2:N, err_mcmc_mean + err_mcmc_std, '--b', 'LineWidth', 1);
plot(2:N, err_ep_mean + err_ep_std, '--r', 'LineWidth', 1);
plot(2:N, err_adf_mean + err_adf_std, '--k', 'LineWidth', 1);  
xlim([2 N])
legend('MCMC','EP','ADF','STDEV');
xlabel('Observations');
ylabel('MI Estimate Error')
yl = ylim;
ylim([0 yl(2)]);

% plot state errs
figure('InvertHardcopy','off','Color',[1 1 1]);
set(gca,'FontSize',14);
hold on
err_mcmc_mean_state = mean(err_mcmc_state, 2);
err_mcmc_std_state = std(err_mcmc_state, 0, 2);
err_random_mean_state = mean(err_random_state, 2);
err_random_std_state = std(err_random_state, 0, 2);
err_ep_mean_state = mean(err_ep_state, 2);
err_ep_std_state = std(err_ep_state, 0, 2);
err_adf_mean_state = mean(err_adf_state, 2);
err_adf_std_state = std(err_adf_state, 0, 2);
plot(1:N, err_random_mean_state, '-g', 'LineWidth', 2);
plot(1:N, err_mcmc_mean_state, '-b', 'LineWidth', 2);
plot(1:N, err_ep_mean_state, '-r', 'LineWidth', 2);
plot(1:N, err_adf_mean_state, '-k', 'LineWidth', 2);  
% plot(1:N, err_random_mean_state - err_random_std_state, '--g', 'LineWidth', 1);
% plot(1:N, err_mcmc_mean_state - err_mcmc_std_state, '--b', 'LineWidth', 1);
% plot(1:N, err_ep_mean_state - err_ep_std_state, '--r', 'LineWidth', 1);
% plot(1:N, err_adf_mean_state - err_adf_std_state, '--k', 'LineWidth', 1);  
plot(1:N, err_random_mean_state + err_random_std_state, '--g', 'LineWidth', 1);
plot(1:N, err_mcmc_mean_state + err_mcmc_std_state, '--b', 'LineWidth', 1);
plot(1:N, err_ep_mean_state + err_ep_std_state, '--r', 'LineWidth', 1);
plot(1:N, err_adf_mean_state + err_adf_std_state, '--k', 'LineWidth', 1);  
xlim([1 N])
legend('Random','MCMC','EP','ADF','STDEV');
xlabel('Observations');
ylabel('State Estimate Error')
yl = ylim;
ylim([0 yl(2)]);
xlim([1, N-1]);


