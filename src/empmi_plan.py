"""
MI planning based on empirical estimate
"""
import numpy as np
import scipy.stats as stats
import utils

def run_empmi_EPLLDA(map_args):
    """
    Wrapper function for empirical MI estimate
    """
    (rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs, words, Nsamp) = map_args
    N = len(docs)

    # compute entropy for each document / word
    MI = np.zeros((N, 1))
    for idx in range(N):
        (MI[idx], rng) = empirical_MI(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs[idx], words[idx], Nsamp)
    return (MI, rng)

def empirical_MI(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, d, n, Nsamp):
    """
    Estimate MI based on empirical mean
    """
    wdn = wordids[d][n]

    # compute cavity parameters of document-topic proportions
    rho_d = gamma[d,:] - zeta[d][n,:]  # array of length K

    # compute cavity for topic parameters
    tau = np.empty_like(lam)
    for k in range(0, K):
        tau[k,:] = lam[k,:] - omega[k][d][n,:]

    # If some components of rho_d or tau_k are non-positive, skip this step
    if (rho_d.min() <= 0 or tau.min() <= 0):
        print("\t\tSkipped (d,n) = (" + str(d) + "," + str(n) \
              + "): min(rho_d) = " + str(rho_d.min()) + ", min(tau) = " + str(tau.min()))

    # Compute mixing proportions
    sum_tau = np.sum(tau, 1)  # array of length K
    pZ = (rho_d * tau[:,wdn]) / sum_tau  # elementwise operation
    pZ = pZ / np.sum( pZ )

    # compute marginal p(y) and conditional distribution p(z|y)
    pZY = pZ[:,np.newaxis] * ppi   # K x Nl
    pY = np.sum( pZY, axis=0 )
    pY = pY / np.sum( pY )  # Nl vector
    pZcondY = pZY / pY[np.newaxis,:]    
    
    # sample  Y / psi from posterior p(.|data)
    y_samp = np.zeros(Nsamp, dtype='int')
    psi_samp = np.zeros((K,W,Nsamp))
    logp_marg_all = np.zeros(Nsamp)
    logp_cond_all = np.zeros(Nsamp)
    for j in range(Nsamp):

        # draw samples
        z = rng.choice(K, p=pZ)
        y = rng.choice(Nl, p=ppi[z,:])

        # sample topics
        psi = np.zeros((K,W))
        logp_psi = np.zeros((K,))
        logp_plusone_psi = np.zeros((K,))        
        for k in range(K):
            lam_k = tau[k,:]
            lam_k_plusone = lam_k.copy()
            lam_k_plusone[wdn] += 1
            if (z==k):
                psi_k = rng.dirichlet( lam_k_plusone )
            else:
                psi_k = rng.dirichlet( lam_k )
            logp_psi[k] = stats.dirichlet.logpdf( psi_k, lam_k )
            logp_plusone_psi[k] = stats.dirichlet.logpdf( psi_k, lam_k_plusone )

        # compute topic marginal p(\psi | data) and condiitonal p(\psi | y, data)
        logp_marg = np.log( pZ )
        logp_cond = np.log( pZcondY )
        for k in range(K):
            this_logp = logp_psi.copy()
            this_logp[k] = logp_plusone_psi[k]
            logp_marg[k] += np.sum( this_logp )
            logp_cond[k] += np.sum( this_logp )
        log_Z_marg = np.max( logp_marg )
        logp_marg_all[j] = log_Z_marg + np.log( np.sum( np.exp( logp_marg - log_Z_marg ) ) )
        log_Z_cond = np.max( logp_cond )
        logp_cond_all[j] = log_Z_cond + np.log( np.sum( np.exp( logp_cond - log_Z_cond ) ) )

    # estimate entropies
    Hmarg = - 1/Nsamp * np.sum( logp_marg_all )
    Hcond = - 1/Nsamp * np.sum( logp_cond_all )

    MI = Hmarg - Hcond
    return (MI, rng)        
