"""
Filtering with Expectation Propagation as primitive (LLDA)
"""

import numpy as np
import scipy.stats as stats
import utils
from scipy.optimize import fmin_l_bfgs_b

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

def run_discvibound_EPLLDA(map_args):
    """
    Wrapper function for "discriminative" variational MI bound
    """
    (rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs, words, Nsamp) = map_args
    N = len(docs)

    # compute entropy for each document / word
    MI = np.zeros((N, 1))
    for idx in range(N):
        (MI[idx], rng) = discvarinf_bound(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs[idx], words[idx], Nsamp)
    return (MI, rng)

def run_estentropy_EPLLDA(map_args):
    """
    Wrapper function for empirical MI estimate
    """
    (rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs, words, Nsamp) = map_args
    N = len(docs)

    # compute entropy for each document / word
    MI = np.zeros((N, 1))
    for idx in range(N):
        (MI[idx], rng) = estimate_entropy_independent(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs[idx], words[idx], Nsamp)
    return (MI, rng)

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

def unwrap(x, K, V):
    w = x[0:K*V]
    w = w.reshape((V,K))
    w0 = x[K*V:]
    return w, w0

def wrap(w, w0, K, V):
    xw = w.reshape(K*V)
    x = np.concatenate((xw, w0))
    return x

def discvar_obj(x, *args):
    """
    Cross-entropy objective for discriminative variational objective
    """
    psi_samp = args[0]
    y_samp = args[1]
    Nsamp = psi_samp.shape[2]
    K = args[2]
    V = args[3]
    w, w0 = unwrap(x, K, V)

    H = 0.
    for i in range(Nsamp):

        # evaluate variational distribution
        psi = psi_samp[...,i] # K x V
        log_q_all = np.diag( np.matmul(psi, w) ) + w0
        log_Z = np.max( log_q_all )
        log_q_all -= log_Z
        log_q_all -= np.log( np.sum( np.exp( log_q_all ) ) )

        # update entropy
        k = y_samp[i]
        H -= 1/Nsamp * log_q_all[k]

    return H

def discvar_grad(x, *args):
    """
    Gradient for cross-entropy objective of discriminative variational objective
    """
    psi_samp = args[0]
    y_samp = args[1]
    Nsamp = psi_samp.shape[2]
    K = args[2]
    V = args[3]
    w, w0 = unwrap(x, K, V)

    # count # of samples for each label Y
    Nsamp_all, junk = np.histogram(y_samp, bins=np.arange(K+1))

    # compute gradient
    gw_p = np.zeros(w.shape) # V x K
    gw_q = np.zeros(w.shape)
    gw0_q = np.zeros(K)
    for i in range(Nsamp):
        k = y_samp[i]
        psi = psi_samp[...,i] # K x V

        # evaluate model gradients for w
        gw_p[:,k] += 1 / Nsamp * psi[k].T

        # evaluate variational distribution
        log_q_all = np.diag( np.dot(psi, w) ) + w0
        log_Z = np.max( log_q_all )
        log_q_all -= log_Z
        log_q_all -= np.log( np.sum( np.exp( log_q_all ) ) )

        # evaluate variational gradients
        gw0_q += 1/Nsamp * np.exp( log_q_all )
        gw_q += 1/Nsamp * np.exp( log_q_all ) * psi.T

    # compute full gradient
    gw = gw_q - gw_p
    gw0 = gw0_q -  Nsamp_all / Nsamp
    g = wrap(gw, gw0, K, V)

    return g


def discvarinf_bound(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, d, n, Nsamp):
    """
    Information theoretic planning based on "discriminative"
    variational bound of mutual information.
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


    # sample Z / Y / psi
    y_samp = np.zeros(Nsamp, dtype='int')
    z_samp = np.zeros(Nsamp, dtype='int')
    psi_samp = np.zeros((K,W,Nsamp))
    for j in range(0, Nsamp):

        # sample discrete components
        z = rng.choice(K, p=C)
        y = rng.choice(Nl, p=ppi[z,:])

        # sample topics
        psi = np.zeros((K,W))
        logp_psi = np.zeros(K)
        logp_plusone_psi = np.zeros(K)
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
            psi[k,:] = psi_k

        # save samples
        z_samp[j] = z
        y_samp[j] = y
        psi_samp[...,j] = psi.copy()

    # optimization
    V = tau.shape[1]
    x, H_cond, d = fmin_l_bfgs_b(discvar_obj, np.zeros(K*V + K), fprime=discvar_grad, args=(psi_samp, y_samp, K, V),
                                  approx_grad=False, iprint=0, factr=1e10)

    # compute marginal entropy
    C_joint = C[:,np.newaxis] * ppi   # K x Nl
    pY = np.sum( C_joint, axis=0 )
    pY = pY / np.sum( pY )  # Nl vector
    H = - np.dot( pY, np.log(pY) )

    MI = H - H_cond
    return (MI, rng)


def estimate_entropy_independent(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, d, n, Nsamp):
    """
    Empirical estimate of conditional entropy H(\psi | Y) by sampling variational posterior
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

    # mixture weights
    sum_tau = np.sum(tau, 1)  # array of length K
    C = (rho_d * tau[:,wdn]) / sum_tau  # elementwise operation
    C = C / np.sum( C )

    # conditional mixture weights
    C_cond = C[:,np.newaxis] * ppi # K x Nl
    C_cond = C_cond / np.sum(C_cond,axis=0)[np.newaxis, :]

    # estimate marginal entropy
    log_pmarg = np.zeros(Nsamp)
    for j in range(0, Nsamp):

        # sample discrete components
        z = rng.choice(K, p=C)
        y = rng.choice(Nl, p=ppi[z,:])

        # sample topics
        psi = np.zeros((K,W))
        logp_psi = np.zeros(K)
        logp_plusone_psi = np.zeros(K)
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
            psi[k,:] = psi_k

        # compute topic conditional probabilities
        logp_vec = np.zeros(K)
        for k in range(K):
            this_logp = logp_psi.copy()
            this_logp[k] = logp_plusone_psi[k]
            logp_vec[k] = np.sum( this_logp )
        log_Z = np.max( logp_vec )
        p_vec = np.exp( logp_vec - log_Z )
        log_pmarg[j] = log_Z + np.log( C.dot(p_vec) )

    # estimate conditional entropy
    log_pcond = np.zeros(Nsamp)
    for j in range(0, Nsamp):

        # sample discrete components
        z = rng.choice(K, p=C)
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
            psi[k,:] = psi_k

        # compute topic conditional probabilities
        logp_vec = np.zeros(K)
        for k in range(K):
            this_logp = logp_psi.copy()
            this_logp[k] = logp_plusone_psi[k]
            logp_vec[k] = np.sum( this_logp )
        log_Z = np.max( logp_vec )
        p_vec = np.exp( logp_vec - log_Z )
        log_pcond[j] = log_Z + np.log( C_cond[:,y].dot(p_vec) )

    # estimate MI
    MI = 1/Nsamp * ( np.sum( log_pcond ) - np.sum( log_pmarg ) )
    return (MI, rng)

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
            psi[k,:] = psi_k

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
    Hcond = 0
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

            # import matplotlib.pyplot as plt
            # fig, axs = plt.subplots(1,2)
            # axs[0].imshow( np.reshape(tau_k, (5,5)), cmap='gray')
            # axs[0].set_title('Cavity')
            # fig.suptitle('wdn: %d, k: %d, l: %d' % (wdn, k,l))
            # axs[1].imshow( np.reshape(E_beta_kl, (5,5)), cmap='gray')
            # axs[1].set_title('Augmented Mean')
            # plt.savefig('/home/pachecoj/Dropbox/Research/journal/042218/tmp/wdn_%d_k_%d_l_%d.png' % (wdn,k,l))
            # plt.close()

            # compute expected entropy
            lambda_kl = (sum(E_beta_kl - E_beta2_kl)) / (sum(E_beta2_kl - E_beta_kl ** 2)) * E_beta_kl
            Hcond += pY[l] * stats.dirichlet.entropy( lambda_kl )

    MI = Hmarg - Hcond
    return (MI, rng)


## OLD CODE ####
#
# def estimate_entropy_dependent(self, d, n, Nsamp):
#     """
#     Empirical estimate of conditional entropy H(\psi | Y)
#     """
#
#     wdn = self._wordids[d][n]
#
#     # compute cavity parameters of document-topic proportions
#     rho_d = self._gamma[d,:] - self._zeta[d][n,:]  # array of length K
#
#     # compute cavity for topic parameters
#     tau = np.empty_like(self._lambda)
#     for k in range(0, self._K):
#         tau[k,:] = self._lambda[k,:] - self._omega[k][d][n,:]
#
#     # If some components of rho_d or tau_k are non-positive, skip this step
#     if (rho_d.min() <= 0 or tau.min() <= 0):
#         print("\t\tSkipped (d,n) = (" + str(d) + "," + str(n) \
#               + "): min(rho_d) = " + str(rho_d.min()) + ", min(tau) = " + str(tau.min()))
#
#     # mixture weights
#     sum_tau = np.sum(tau, 1)  # array of length K
#     C = (rho_d * tau[:,wdn]) / sum_tau  # elementwise operation
#     C = C / np.sum( C )
#
#     # conditional mixture weights
#     C_cond = C[:,np.newaxis] * self._ppi # K x Nl
#     C_cond = C_cond / np.sum(C_cond,axis=0)[np.newaxis, :]
#
#     # import ipdb;  ipdb.set_trace()
#
#     log_pcond = np.zeros(Nsamp)
#     log_pmarg = np.zeros(Nsamp)
#     for j in range(0, Nsamp):
#
#         # sample discrete components
#         z = np.random.choice(self._K, p=C)
#         y = np.random.choice(self._Nl, p=self._ppi[z,:])
#
#         # sample topics
#         psi = np.zeros((self._K,self._W))
#         logp_psi = np.zeros((self._K,))
#         logp_plusone_psi = np.zeros((self._K,))
#         for k in range(self._K):
#             lam_k = tau[k,:]
#             lam_k_plusone = lam_k.copy()
#             lam_k_plusone[wdn] += 1
#             if (z==k):
#                 psi_k = np.random.dirichlet( lam_k_plusone )
#             else:
#                 psi_k = np.random.dirichlet( lam_k )
#             logp_psi[k] = stats.dirichlet.logpdf( psi_k, lam_k )
#             logp_plusone_psi[k] = stats.dirichlet.logpdf( psi_k, lam_k_plusone )
#             psi[k,:] = psi_k
#
#         # compute topic conditional probabilities
#         logp_vec = np.zeros(self._K)
#         for k in range(self._K):
#             this_logp = logp_psi.copy()
#             this_logp[k] = logp_plusone_psi[k]
#             logp_vec[k] = np.sum( this_logp )
#         log_Z = np.max( logp_vec )
#         p_vec = np.exp( logp_vec - log_Z )
#         log_pcond[j] = log_Z + np.log( C_cond[:,y].dot(p_vec) )
#         log_pmarg[j] = log_Z + np.log( C.dot(p_vec) )
#
#     # estimate MI
#     MI = 1/Nsamp * ( np.sum( log_pcond ) - np.sum( log_pmarg ) )
#     return MI
