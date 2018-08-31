%% Compile Drovandi et al. 2014 code

% compile the c code for MH step for exponential model (this only needs to be done once)
mex CFLAGS='-std=c99 -fPIC' posterior_dist_rw_smc_exp_mex.c mt19937ar.c;

% compile the c code for MH step for power model (this only needs to be done once)
mex CFLAGS='-std=c99 -fPIC' posterior_dist_rw_smc_power_mex.c mt19937ar.c;

% compile the code for performing systematic resampling (this only needs to be done once)
mex CFLAGS='-std=c99 -fPIC' systematic_resampling.c mt19937ar.c;
