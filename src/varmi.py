"""
Filtering with Expectation Propagation as primitive (LLDA)
"""

import numpy as np
import scipy.stats as stats
import utils

def run_vibound_EPLLDA(map_args):
    """
    Wrapper function for variational MI bound
    """
    (W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs, words) = map_args
    N = len(docs)

    # compute entropy for each document / word
    MI = np.zeros((N, 1))
    for idx in range(N):
        MI[idx] = varinf_bound(W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs[idx], words[idx])
    return MI

def run_estMIhybrid_EPLLDA(map_args):
    """
    Wrapper function for hybrid MI estimate
    """
    (rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs, words, Nsamp) = map_args
    N = len(docs)

    # compute entropy for each document / word
    MI = np.zeros((N, 1))
    for idx in range(N):
        (MI[idx], rng) = varinf_hybrid_bound(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs[idx], words[idx], Nsamp)
    return (MI, rng)

def varinf_bound(W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, d, n):
    """
    Information theoretic planning based on variational bound of
    mutual information.
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
    C = (rho_d * tau[:,wdn]) / sum_tau  # elementwise operation
    C = C / np.sum( C )

    C_cond = C[:,np.newaxis] * ppi   # K x Nl
    pY = np.sum( C_cond, axis=0 )
    pY = pY / np.sum( pY )  # Nl vector
    C_cond = C_cond / pY[np.newaxis,:]

    # store marginal entropy of topic components
    H_psi = np.zeros(K)
    H_psi_plusone = np.zeros(K)
    for k in range(K):
        lam_k = tau[k,:]
        lam_k_plusone = lam_k.copy()
        lam_k_plusone[wdn] += 1
        H_psi[k] = stats.dirichlet.entropy(lam_k)
        H_psi_plusone[k] = stats.dirichlet.entropy(lam_k_plusone)

    # compute MI bound
    Hcond = 0
    Hmarg = 0
    for k in range(K):
        tau_k = tau[k,:]
        sum_tau_k = sum_tau[k]
        tau_k_wdn = tau_k[wdn]
        C_k = C[k]

        # marginal entropy
        this_H_vec = H_psi.copy()
        this_H_vec[k] = H_psi_plusone[k]
        Hmarg += C_k * sum( this_H_vec )

        # conditional entropy
        for l in range(Nl):
            C_kl = C_cond[k,l]

            # compute conditional moments
            E_beta_kl = \
              tau_k / (sum_tau_k + 1) + (1 - C_kl) * tau_k / (sum_tau_k * (sum_tau_k + 1))
            E_beta_kl[wdn] += C_kl / (sum_tau_k + 1)
            E_beta2_kl = tau_k * (tau_k + 1) / ((sum_tau_k + 1) * (sum_tau_k + 2)) \
              + 2 * (1 - C_kl) * tau_k * (tau_k + 1) / \
              (sum_tau_k * (sum_tau_k + 1) * (sum_tau_k + 2))
            E_beta2_kl[wdn] += 2 * C_kl * (tau_k_wdn + 1) / \
              ((sum_tau_k + 1) * (sum_tau_k + 2))

            # compute expected entropy
            lambda_kl = (sum(E_beta_kl - E_beta2_kl)) / (sum(E_beta2_kl - E_beta_kl ** 2)) * E_beta_kl
            Hcond += pY[l] * stats.dirichlet.entropy( lambda_kl )

    return Hmarg - Hcond

def varinf_hybrid_bound(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, d, n, Nsamp):
    """
    Information theoretic planning based on variational bound of
    mutual information.  This version is a hybrid which estimates the marginal
    entropy H(\psi) through samples of the variational posterior, and bounds the
    conditional entropy H(\psi | Y)
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
    C = (rho_d * tau[:,wdn]) / sum_tau  # elementwise operation
    C = C / np.sum( C )

    # sample-based estimate of marginal entropy
    log_pmarg = np.zeros(Nsamp)
    for j in range(0, Nsamp):

        # sample topic
        z = rng.choice(K, p=C)
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

        # compute topic conditional probabilities
        logp_vec = np.zeros(K)
        for k in range(K):
            this_logp = logp_psi.copy()
            this_logp[k] = logp_plusone_psi[k]
            logp_vec[k] = np.sum( this_logp )
        log_Z = np.max( logp_vec )
        p_vec = np.exp( logp_vec - log_Z )
        log_pmarg[j] = log_Z + np.log( C.dot(p_vec) )

    # estimate marginal entropy
    Hmarg = - 1/Nsamp * np.sum( log_pmarg )

    # conditional mixture weights
    C_cond = C[:,np.newaxis] * ppi   # K x Nl
    pY = np.sum( C_cond, axis=0 )
    pY = pY / np.sum( pY )  # Nl vector
    C_cond = C_cond / pY[np.newaxis,:]

    # bound conditional entropy
    Hcond = 0.
    for k in range(0, K):
        tau_k = tau[k,:]
        sum_tau_k = sum_tau[k]
        tau_k_wdn = tau_k[wdn]
        for l in range(0, Nl):
            C_kl = C_cond[k,l]

            # compute conditional moments
            E_beta_kl = \
              tau_k / (sum_tau_k + 1) + (1 - C_kl) * tau_k / (sum_tau_k * (sum_tau_k + 1))
            E_beta_kl[wdn] += C_kl / (sum_tau_k + 1)
            E_beta2_kl = tau_k * (tau_k + 1) / ((sum_tau_k + 1) * (sum_tau_k + 2)) \
              + 2 * (1 - C_kl) * tau_k * (tau_k + 1) / \
              (sum_tau_k * (sum_tau_k + 1) * (sum_tau_k + 2))
            E_beta2_kl[wdn] += 2 * C_kl * (tau_k_wdn + 1) / \
              ((sum_tau_k + 1) * (sum_tau_k + 2))

            # compute expected entropy
            lambda_kl = (sum(E_beta_kl - E_beta2_kl)) / (sum(E_beta2_kl - E_beta_kl ** 2)) * E_beta_kl
            Hcond += pY[l] * stats.dirichlet.entropy( lambda_kl )

    MI = Hmarg - Hcond
    return (MI, rng)
