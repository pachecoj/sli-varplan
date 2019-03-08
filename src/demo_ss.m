% DEMO_SS - Demonstrate sensor selection for the nonlinear state-space model in
%   Example 1 of:
%
%   O. Cappe, S. Godsill, & E. Moulines, "An Overview of Existing Methods 
%     and Recent Advances in Sequential Monte Carlo". Proceedings of the 
%     IEEE, vol. 95, pp. 899-924, 2007%
%
% J. Pacheco, 2013, 2017
%

%% NOTES:20
% * Rename existing "optimal" to "oracle"
% * "Optimal" Sensor selection using numerical approximation of mutual inf.

clear variables
close all

seed = 1;
rng(seed);

% method defs
ORACLE_OPT = 1;
OPT_OPT = 2;
VI_OPT = 3;
ORACLE_PF = 4;
MCMC_PF = 5;
VI_PF = 6;
NMETHODS = 6;

% model params
sig_u = 10;       % process noise
alpha = 0.05;      % measurement variance factor
T = 20;           % # scans
mu_0 = 0;         % prior mean
sig_0 = sig_u;    % prior stdev
N = 20;            % number of random trials
Ngrid = 1000;     % # discrete grid points for numerical approx.
numSensors = 10;      % # of sensor locations
num_particles = 500;  % # of particle filter particles

% generate state & place sensors
lb = -100;
ub = 100;
loc = linspace(lb, ub, numSensors)';
x_grid = linspace(lb, ub, Ngrid);
dx = x_grid(2) - x_grid(1);

% sample state-space model
rms_state = NaN(NMETHODS,N);
rms_sensor = NaN(NMETHODS,N);
for n=1:N  
  fprintf('Run %d...', n);

  % simulate data / evaluate model
  [x, y_all] = simtrack( T, loc, alpha, mu_0, sig_0, sig_u );    
  [prior, trans, like_all] = eval_model(x_grid, loc, y_all, mu_0, sig_0, sig_u, alpha);
  
  % oracle sensor selection
  D = abs( bsxfun(@minus, loc, x') );
  [~,Ks_oracle] = min(D, [], 1);
  
  % segment oracle likelihoods / observations
  like_oracle = cell(T,1);
  y_oracle = NaN(T,1);
  for t=1:T
    like_oracle{t} = like_all{t}(:,Ks_oracle(t));
    y_oracle(t) = y_all(t,Ks_oracle(t));
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SELECTION: ORACLE, INFERENCE: OPTIMAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Oracle selection means the sensor closest to the actual target
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  log_marg = sum_product_hmm( prior, trans, like_oracle );
  [~,idx] = max( log_marg, [], 1 );  
  mpm_oracle = x_grid(idx);    
%   rms_state(VI_OPT,n) = sqrt(mean((mpm_hmmvi - mpm_opt).^2));
    
  % show true marginal
  if N==1
    figure('InvertHardcopy','off','Color',[1 1 1]);
    subplot(2,3,1);
    [POS,TIME] = meshgrid(x_grid, (1:(T+1))-0.5);
    s = pcolor( TIME, POS, [-exp(log_marg), zeros(numel(x_grid),1)]' );
    hold on;
    colormap gray
    set(s,'EdgeColor','none');
    set(gca,'FontSize',14);
    plot( 1:T, x, '-og');
    plot(1:T, mpm_oracle, '--or');
    scatter(1:T, loc( Ks_oracle ), 'ok')
    legend('Marginal','True State','MPM','Selection');
    xlabel('Time');
    ylabel('Position');
    title('Selection: Oracle, Inference: Exact');
    ylim([-50 50]);
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SELECTION: OPTIMAL, INFERENCE: OPTIMAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Optimal selection is the sensor maximizing MI by
  % numerical integration.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [log_marg, Ks_opt] = sum_product_hmm_optselect( ...
    x_grid, prior, trans, like_all, loc, alpha, 1 );
  [~,idx] = max( log_marg, [], 1 );
  mpm_opt = x_grid(idx);
  
  % show true marginal  
  if N==1
    subplot(2,3,2);
    [POS,TIME] = meshgrid(x_grid, (1:(T+1))-0.5);
    s = pcolor( TIME, POS, [-exp(log_marg), zeros(numel(x_grid),1)]' );
    hold on;
    colormap gray
    set(s,'EdgeColor','none');
    set(gca,'FontSize',14);
    plot( 1:T, x, '-og');
    plot(1:T, mpm_opt, '--or');
    scatter(1:T, loc( Ks_opt ), 'ok')
    xlabel('Time');
    ylabel('Position');
    title('Selection: Exact MI, Inference: Exact');
    ylim([-50 50]);
  end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SELECTION: VARIATIONAL, INFERENCE: OPTIMAL
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [log_marg, Ks_hmmvi] = sum_product_hmm_viselect( ...
    x_grid, prior, trans, like_all, loc, alpha, 1 );
  [~,idx] = max( log_marg, [], 1 );
  mpm_hmmvi = x_grid(idx);
  rms_state(VI_OPT,n) = sqrt(mean((mpm_hmmvi - mpm_opt).^2));
  rms_sensor(VI_OPT,n) = sqrt(mean((loc(Ks_opt) - loc(Ks_hmmvi)).^2));
  
  % plot optimal marginals w/ VI selection
  if N==1
    subplot(2,3,3);
    [POS,TIME] = meshgrid(x_grid, (1:(T+1))-0.5);
    s = pcolor( TIME, POS, [-exp(log_marg), zeros(numel(x_grid),1)]' );
    hold on;
    colormap gray
    set(s,'EdgeColor','none');
    set(gca,'FontSize',14);
    plot( 1:T, x, '-og');
    plot(1:T, mpm_hmmvi, '--or');
    scatter(1:T, loc( Ks_hmmvi ), 'ok')
    xlabel('Time');
    ylabel('Position');
    title('Selection: VIP, Inference: Exact');
    ylim([-50 50]);
  end
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SELECTION: ORACLE, INFERENCE: PARTICLE FILTER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  kdsig = 1;
  resamp_int = 1;  
  [x_pf,w_pf] = particle_filter( ...
    num_particles, resamp_int, Ks_oracle, loc, y_oracle, mu_0, sig_0, sig_u, alpha  );
  marg_pf = kernel_density_gauss( x_grid, x_pf, w_pf, kdsig );
  [~,idx] = max( marg_pf, [], 1 );
  mpm_pf = x_grid(idx);
      
  % plot PF marginals w/ oracle selection
  if N==1
    subplot(2,3,4);
    s = pcolor( TIME, POS, [-marg_pf, zeros(numel(x_grid),1)]' );
    colormap gray
    set(s,'EdgeColor','none');
    set(gca,'FontSize',14);
    hold on;
    t_mat = repmat(1:T, [num_particles, 1]);
    scatter(t_mat(:), x_pf(:), 1, '.k');
    plot( 1:T, x, '-og');
    plot(1:T, mpm_pf, '--or');
    scatter(1:T, loc( Ks_oracle ), 'ok')
    legend('','Particles');
    xlabel('Time');
    ylabel('Position');
    title('Selection: Oracle, Inference: PF');
    ylim([-50 50]);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SELECTION: MCMI, INFERENCE: PARTICLE FILTER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [x_pfmc, Ks_pfmc, y_pfmc, w_pfmc] = pf_mcsel( ...
    num_particles, resamp_int, loc, y_all, mu_0, sig_0, sig_u, alpha  );
  marg_pfmc = kernel_density_gauss( x_grid, x_pfmc, w_pfmc, kdsig );
  [~,idx] = max( marg_pfmc, [], 1 );
  mpm_pfmc = x_grid(idx);
  rms_state(MCMC_PF,n) = sqrt(mean((mpm_pfmc - mpm_opt).^2));
  rms_sensor(MCMC_PF,n) = sqrt(mean((loc(Ks_opt) - loc(Ks_pfmc)).^2));
  
  % plot PF marginals w/ MCMC selection
  if N==1
    subplot(2,3,5);
    s = pcolor( TIME, POS, [-marg_pfmc, zeros(numel(x_grid),1)]' );
    colormap gray
    set(s,'EdgeColor','none');
    set(gca,'FontSize',14);
    hold on;
    plot( 1:T, x, '-og');
    t_mat = repmat(1:T, [num_particles, 1]);
    scatter(t_mat(:), x_pfmc(:), 1, '.k');
    plot(1:T, mpm_pfmc, '--or');
    scatter(1:T, loc( Ks_pfmc ), 'ok')
    xlabel('Time');
    ylabel('Position');
    title('Selection: Empirical, Inference PF');
    ylim([-50 50]);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% SELECTION: VARIATIONAL, INFERENCE: PARTICLE FILTER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [x_pfvi, Ks_pfvi, y_pfvi, w_pfvi] = pf_varinf( ...
    num_particles, resamp_int, loc, y_all, mu_0, sig_0, sig_u, alpha  );
  marg_pfvi = kernel_density_gauss( x_grid, x_pfvi, w_pfvi, kdsig );
  [~,idx] = max( marg_pfvi, [], 1 );
  mpm_pfvi = x_grid(idx);
  rms_state(VI_PF,n) = sqrt(mean((mpm_pfvi - mpm_opt).^2));
  rms_sensor(VI_PF,n) = sqrt(mean((loc(Ks_opt) - loc(Ks_pfvi)).^2));
  
  % plot PF marginals w/ VI selection
  if N==1
    subplot(2,3,6);
    s = pcolor( TIME, POS, [-marg_pfvi, zeros(numel(x_grid),1)]' );
    colormap gray
    set(s,'EdgeColor','none');
    set(gca,'FontSize',14);
    hold on;
    plot( 1:T, x, '-og');
    t_mat = repmat(1:T, [num_particles, 1]);
    scatter(t_mat(:), x_pfvi(:), 1, '.k');
    plot(1:T, mpm_pfvi, '--or');
    scatter(1:T, loc( Ks_pfvi ), 'ok')
    xlabel('Time');
    ylabel('Position');
    title('Selection: VIP, Inference PF');
    ylim([-50 50]);
  end
  
  fprintf('done.\n');
end

% ORACLE_OPT = 1;
% OPT_OPT = 2;
% VI_OPT = 3;
% ORACLE_PF = 4;
% MCMC_PF = 5;
% VI_PF = 6;

% plot MPM error
figure('InvertHardcopy','off','Color',[1 1 1]);
set(gca,'FontSize',14);
boxplot(rms_state([ORACLE_OPT, OPT_OPT, VI_OPT],:)');
set(gca,'XTickLabel',{'Oracle / Optimal','PF VI','MCMI'});
ylabel('RMS Error');
title('State Error (Exact Inference)');
figure('InvertHardcopy','off','Color',[1 1 1]);
set(gca,'FontSize',14);
boxplot(rms_state([ORACLE_PF, MCMC_PF, VI_PF],:)');
set(gca,'XTickLabel',{'Oracle / Optimal','PF VI','MCMI'});
ylabel('RMS Error');
title('State Error (Particle Filter)');

% plot selection error
figure('InvertHardcopy','off','Color',[1 1 1]);
set(gca,'FontSize',14);
boxplot(rms_sensor([ORACLE_OPT, OPT_OPT, VI_OPT, ORACLE_PF, MCMC_PF, VI_PF],:)');
set(gca,'XTickLabel',{'Exact VI','PF VI','MCMI'});
ylabel('RMS Error');
title('RMS Sensor Selection Error');
