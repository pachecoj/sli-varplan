"""
Discriminative variational MI planning
"""

import numpy as np
import scipy.stats as stats
import utils
import time
from scipy.optimize import fmin_l_bfgs_b

def run_discvibound_EPLLDA(map_args):
    """
    Wrapper function for "discriminative" variational MI bound
    """
    (rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs, words, Nsamp) = map_args
    N = len(docs)

    # compute entropy for each document / word
    MI = np.zeros(N)
    t = np.zeros(N)
    for idx in range(N):
        (MI[idx], rng) = discvarinf_bound(rng, W, K, Nl, ppi, wordids, gamma, zeta, lam, omega, docs[idx], words[idx], Nsamp)
    return (MI, rng)

def unwrap(x, K, V):
    w = x[0:K*K]
    w = w.reshape((K,K))
    w0 = x[K*K:]
    return w, w0

def wrap(w, w0, K, V):
    xw = w.reshape(K*K)
    x = np.concatenate((xw, w0))
    return x

def discvar_obj(x, *args):
    """
    Cross-entropy objective for discriminative variational objective
    """
    psi_samp = args[0]
    y_samp = args[1]
    Nsamp = psi_samp.shape[1]
    K = args[2]
    V = args[3]
    w, w0 = unwrap(x, K, V)

    # Parameter w is KxK where w[:,i] is a K-vector
    # for i^th label p(Y=i | psi)

    H = 0.
    for i in range(Nsamp):

        # evaluate variational distribution
        psi = psi_samp[...,i] # K-vector
        log_q_all = np.dot(w.T,psi) + w0
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
    Nsamp = psi_samp.shape[1]
    K = args[2]
    V = args[3]
    w, w0 = unwrap(x, K, V)

    # count # of samples for each label Y
    Nsamp_all, junk = np.histogram(y_samp, bins=np.arange(K+1))

    # compute gradient
    gw_p = np.zeros(w.shape) # K x K
    gw_q = np.zeros(w.shape)
    gw0_q = np.zeros(K)
    for i in range(Nsamp):
        k = y_samp[i]
        psi = psi_samp[...,i] # K vector

        # evaluate model gradients for w
        gw_p[:,k] += 1 / Nsamp * psi

        # evaluate variational distribution
        log_q_all = np.dot(w.T,psi) + w0
        log_Z = np.max( log_q_all )
        log_q_all -= log_Z
        log_q_all -= np.log( np.sum( np.exp( log_q_all ) ) )

        # evaluate variational gradients
        gw0_q += 1/Nsamp * np.exp( log_q_all )
        gw_q += 1/Nsamp * np.exp( log_q_all ) * psi

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
    pZ = (rho_d * tau[:,wdn]) / sum_tau  # elementwise operation
    pZ = pZ / np.sum( pZ )

    # sample Z / Y / psi
    y_samp = np.zeros(Nsamp, dtype='int')
    z_samp = np.zeros(Nsamp, dtype='int')
    psi_samp = np.zeros((K,W,Nsamp))
    for j in range(0, Nsamp):

        # sample discrete components
        z = rng.choice(K, p=pZ)
        y = rng.choice(Nl, p=ppi[z,:])

        # sample topics
        psi = np.zeros((K,W))
        for k in range(K):
            lam_k = tau[k,:]
            lam_k_plusone = lam_k.copy()
            lam_k_plusone[wdn] += 1
            if (z==k):
                psi_k = rng.dirichlet( lam_k_plusone )
            else:
                psi_k = rng.dirichlet( lam_k )
            psi[k,:] = psi_k

        # save samples
        z_samp[j] = z
        y_samp[j] = y
        psi_samp[...,j] = psi.copy()

    # isolate wdn^th word from topics
    psi_samp = psi_samp[:,wdn,:] # K x Nsamp

        
    # optimization
    V = tau.shape[1]
    x, H_cond, d = fmin_l_bfgs_b(discvar_obj, np.zeros(K*K + K), fprime=discvar_grad, args=(psi_samp, y_samp, K, V),
                                  approx_grad=False, iprint=1, factr=1e10)

    # compute marginal entropy
    pZY = pZ[:,np.newaxis] * ppi   # K x Nl
    pY = np.sum( pZY, axis=0 )
    pY = pY / np.sum( pY )  # Nl vector
    H_marg = - np.dot( pY, np.log(pY) )

    import IPython
    IPython.embed()

    MI = H_marg - H_cond
    return (MI, rng)
