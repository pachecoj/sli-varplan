function x_samp = mcmc_gmm(Nsamp, model, x, this_y, k_vec, x_grid, fnplot, h_fig)
% MCMC_GMM - MCMC sampler for GMM
%
% J. Pacheco, 2017

  % CHEAT : Build dense approximation of true posterior
  [true_posterior_pdf, ~, best_fit_pdf] = true_posterior_gmm( model, this_y, x_grid, k_vec );
  
  % draw samples
  dx = x_grid(2) - x_grid(1);
  cdf = cumsum(true_posterior_pdf) * dx;
  [~, idx_samp] = histc( rand(Nsamp,1), cdf);
  x_samp = x_grid( idx_samp );
  
end