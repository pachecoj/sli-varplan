"""
Filtering with Gibbs sampler as primitive (LLDA)
"""

import copy, utils, math, time
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def run_estEntropyGibbs(map_args):

    rng, phi, z, count_z, ppi, K, Nl, Nd, beta, docs, words = map_args
    N = len(docs)

    # estimate entropies
    Hhat = np.zeros((N, 1))
    for idx in range(N):
        d = docs[idx]
        n = words[idx]
        this_z = z[n,d,:]
        (Hhat[idx], rng) = estimate_entropy_single(rng, phi, this_z, count_z, ppi, K, Nl, Nd, beta, d, n)
    return (Hhat, rng)

def init_samples(rng, W, K, D, alpha, beta):
    phi = np.zeros((K,W))
    for k in range(K):
        phi[k,:] = rng.dirichlet(beta[k,:])
    theta = rng.dirichlet(alpha * np.ones(K), size=D )
    return (phi, theta)

def gibbs(map):
    """
    Gibbs sampler for LLDA
    """

    (rng, W, K, D, Nd, Nl, labels, wordids, alpha, beta, ppi, Nsamp, burn) = map

    # samples
    phi_all = np.zeros((K,W,Nsamp))
    theta_all = np.zeros((D,K,Nsamp))
    z_all = np.zeros((Nd,D,Nsamp), dtype='int')
    count_z_all = np.zeros((W,K,Nsamp), dtype='int')

    logp_map = - np.Inf
    phi_map = np.zeros((K,W))
    logp_trace = np.zeros((Nsamp,burn))
    for isamp in range(Nsamp):
        phi, theta = init_samples(rng, W, K, D, alpha, beta)
        z_label = np.zeros((Nd,D), dtype='int')

        # sample variables tot_iter times
        for iter in range(burn):

            # sample topic labels (Z)
            count_z = np.zeros((W, K), dtype='int')
            for d in range(D):
                n_train = np.where( labels[d,:] > -1 )[0]
                wdn_all = wordids[d]

                # topic assignment distribution
                logp = np.log( theta[d,:][:,np.newaxis] ) + np.log( phi[:,wdn_all] ) # K x Nd
                if n_train.size:
                    this_ppi = ppi[labels[d,n_train],:]
                    logp[:,n_train] += np.log( this_ppi.T )
                Z = np.max(logp, axis=0)
                p = np.exp( logp - Z[np.newaxis,:] )
                p = p / np.sum(p, axis=0)[np.newaxis,:]

                # sample topic assignments (SLOW)
                z_d = np.zeros((Nd,K), dtype='int')
                for n in range(Nd):
                    z_d[n,:] = rng.multinomial(1, p[:,n])

                # sample document-level proportions (\theta)
                theta[d,:] = rng.dirichlet( alpha + np.sum(z_d, axis=0) )

                # running count of assignments (SLOW)
                for n in range(Nd):
                    count_z[wdn_all[n],:] += z_d[n,:]

                # convert Z's to label matrix (Nd x D)
                z_label[:,d] = np.sum( z_d * np.arange(K), axis=1)

            # sample topics (\phi)
            for k in range(K):
                phi[k,:] = rng.dirichlet( beta[k,:] + count_z[:,k] )

            # eval log-probability (save MAP estimate)
            logp_new  = eval_logjoint(alpha, beta, ppi, labels, wordids, phi, theta, z=z_label.T)
            if logp_new > logp_map:
                phi_map = phi.copy()
                logp_map = logp_new.copy()
            logp_trace[isamp,iter] = logp_new

        # save samples
        phi_all[:,:,isamp] = phi
        theta_all[:,:,isamp] = theta
        z_all[:,:,isamp] = z_label
        count_z_all[:,:,isamp] = count_z

    # return samples
    return (phi_map, logp_map, phi_all, theta_all, z_all, count_z_all, logp_trace, rng)

def eval_logjoint(alpha, beta, ppi, labels, wordids, phi, theta, z=[]):
    """
    Evaluates log-pdf of the joint for topics phi and topic proportions theta.
    Marginalizes out topic assignments.
    """

    K, W = phi.shape
    D = len( wordids )
    Nd = len( wordids[0] )

    # import ipdb;  ipdb.set_trace()

    # topic prior
    logp = 0.
    for k in range(K):
        logp += stats.dirichlet.logpdf(phi[k,:], beta[k,:])

    # documents
    for d in range(D):

        # topic proportion prior
        logp += stats.dirichlet.logpdf(theta[d,:], alpha * np.ones(K))

        # words
        for n in range(Nd):
            wdn = wordids[d][n]

            # marginalize over topic assignments
            if len(z)==0:
                logp_z = np.zeros(K)
                for k in range(K):
                    logp_z[k] = np.log( theta[d,k] ) + np.log( phi[k,wdn] )
                    l = labels[d,n]
                    if l > -1:
                        logp_z[k] += np.log( ppi[k,l] )

                # numerically stable marginalization
                log_Z = np.max(logp_z)
                logp += log_Z + np.log( np.sum( np.exp( logp_z - log_Z ) ) )

            # compute topic assignment probabilities
            else:
                logp += np.log( theta[d,z[d,n]] ) + np.log( phi[z[d,n],wdn] )
                l = labels[d,n]
                if l > -1:
                    logp += np.log( ppi[z[d,n],l] )

    return logp

def estimate_entropy_single(rng, phi, z, count_z, ppi, K, Nl, Nd, beta, d, n):

    # split samples
    Nsamp_tot = len(z)
    Nsamp = int(np.ceil(Nsamp_tot / 2))
    z_samp = z[0:Nsamp]
    z_samp_marg = z[Nsamp:]
    count_z_samp_marg = count_z[:, :, Nsamp:]  # W x K x Nsamp
    phi_samp = phi[:, :, 0:Nsamp] # K, W, Nsamp

    # sample Y's from joint
    likelihood = ppi[z_samp, :]  # Nsamp x Nl
    y_samp = np.zeros((Nsamp, Nl), dtype='int')
    for i in range(len(y_samp)): # TODO: vectorize
        y_samp[i, :] = rng.multinomial(1, likelihood[i, :])

    # estimate marginal posterior p(y | Data)
    likelihood_marg = ppi[z_samp_marg, :]  # Nsamp x Nl
    phat_y = np.zeros(Nsamp)
    for i in range(Nsamp): # TODO: vectorize
        this_y_samp = y_samp[i, :]  # Nl-vec
        this_likelihood_marg = likelihood_marg * this_y_samp[np.newaxis, :]
        this_likelihood_marg = np.sum(this_likelihood_marg, axis=1)  # Nsamp-vec
        phat_y[i] = 1 / Nsamp * np.sum(this_likelihood_marg, axis=0)

    # compute Dirichlet conditional p(\phi | z, words)
    log_pcond_phi = np.zeros((Nsamp, Nsamp))  # 1st dim is phi-samples / 2nd dim is z-samples
    beta_sim = beta.T[:, :, np.newaxis] + count_z_samp_marg  # W x K x Nsamp
    for k in range(K):
        log_pcond_phi += \
            utils.dirichlet_logpdf(phi_samp[k, ...], beta_sim[:, k, :])

    # estimate joint p(\phi, y | words)
    log_phat_joint = np.zeros(Nsamp)
    for i in range(Nsamp):  # TODO: vectorize
        this_y_samp = y_samp[i, :]  # Nl-vec
        this_likelihood_marg = likelihood_marg * this_y_samp[np.newaxis, :] # Nsamp (z's) x Nl
        this_likelihood_marg = np.sum(this_likelihood_marg, axis=1)  # Nsamp-vec (z's)
        tmp_logp = np.log(this_likelihood_marg) + log_pcond_phi[i, :] # Nsamp-vec (z's)

        # numerically stable average
        log_Z = np.max(tmp_logp)
        tmp_logp = tmp_logp - log_Z
        log_phat_joint[i] = log_Z + np.log(np.sum(np.exp(tmp_logp), axis=0))

    # estimate conditional entropy
    Hhat = 1 / Nsamp * np.sum(np.log(phat_y) - log_phat_joint, axis=0)
    return (Hhat, rng)

