# batch.py: Code for running standard (offline) batch variational inference

# This code suite is largely adapted from the online VB (aka stochastic
# variational Bayes) code of
# Matthew D. Hoffman, Copyright (C) 2010
# found here: http://www.cs.princeton.edu/~blei/downloads/onlineldavb.tar
# and also of
# Chong Wang, Copyright (C) 2011
# found here: http://www.cs.cmu.edu/~chongw/software/onlinehdp.tar.gz
#
# Adapted by: Nick Boyd, Tamara Broderick, Andre Wibisono, Ashia C. Wilson
#
# This program is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, eith
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details. You should have received a copy of the GNU General
# Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.

# import matplotlib.pyplot as plt
# from multiprocessing import Pool
import pickle, string, getopt, sys, random, time, re, pprint, copy, math
import utils, gibbs_llda, ep_llda
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

# Filtering with EP as primitive (LLDA)
class LLDA_EP:

    def __init__(self, pool, rngs, W, K, alpha, beta, ppi, Nl, labels, wordids, maxiters, \
                 threshold, numWorkers, useNewton, Nsamp, silent=False):
        """
        Arguments:
        K: Number of topics
        W: Total number of words in the vocab
        alpha: Hyperparameter for prior on weight vectors theta
        beta: Hyperparameter for prior on topics beta
        Nl: Number of labels
        labels: Initial word labels
        wordids: list of lists with word indices
        maxiters: Max number of iterations to allow to converge
        threshold: Threshold for convergence
        numWorkers: Number of concurrent workers for parallel planning
        useNewton: Boolean, if true we use Newton's method for solving Dirichlet moment-matching, else use approximate MLE
        silent: (optional) Silence verbose output
        """

        self._str = "LLDA_EP_%d_%g_%r" % (maxiters,threshold,useNewton)
        self._K = K
        self._W = W
        self._alpha = alpha
        self._ppi = ppi
        self._Nl = Nl
        self._labels = labels
        self._maxiters = maxiters
        self._threshold = threshold
        self._numWorkers = numWorkers
        self._useNewton = useNewton
        self._silent = silent
        self._D = len(wordids)
        self._wordids = wordids
        self._zeta = list()
        self._omega = list()
        self._lambda = np.empty(0)
        self._gamma = np.empty(0)
        if np.isscalar(beta):
            self._beta = copy.deepcopy(beta) * np.ones((self._K, self._W))
        else:
            self._beta = copy.deepcopy(beta)

        # Get length of each document
        self._N = [len(ids) for ids in self._wordids]
        self._Nsamp = Nsamp

        # random worker pool
        self._pool = pool
        self._rngs = rngs


    def init_messages(self):

        # Initialize the variational distribution q(theta | gamma)
        # with gamma_d = alpha + \sum_n zeta_dn
        # by setting all zeta_dn = 0.
        # zeta is a list of D matrices, with zeta[d] being an N_d-by-K matrix,
        # where N_d is the number of words in the d-th document
        self._zeta.clear()
        for d in range(0, self._D):
            self._zeta.append(np.zeros([self._N[d], self._K]))

        # Initialize the variational distribution q(beta | lambda)
        # with lambda_k = eta_k + \sum_d \sum_n omega_kdn
        # by setting all omega_kdn = 0.
        # omega is a list of K entries, with omega[k] being a list of D matrices,
        # with omega[k][d] being an N_d-by-W matrix.
        self._omega.clear()
        for k in range(0, self._K):
            omega_k = list()
            for d in range(0, self._D):
                omega_k.append(np.zeros([self._N[d], self._W]))
            self._omega.append(omega_k)

        # Initialize lambda and gamma to be equal to the prior beta and alpha
        # (Note: This assumes omega and zeta to be initialized to 0)
        self._lambda = self._beta.copy()
        self._gamma = self._alpha * np.ones([self._D, self._K])

    def label_element(self, labels_all, algname):
        """
        Select word to label
        """
        doclist, wordlist = np.where( self._labels == -1 )
        sizeOfChunks = int(math.ceil(len(doclist) / float(self._numWorkers)))
        docchunks = utils.chunk(doclist, sizeOfChunks)
        wordchunks = utils.chunk(wordlist, sizeOfChunks)
        numTasks = len(docchunks)

        # random planning
        if (algname == 'random'):
            N = len( doclist )
            idx = np.random.randint( N )
            d_new, n_new = (doclist[idx], wordlist[idx])
            MI = []

        # MI-based planning
        else:

            # full variational
            if (algname == 'variational'):
                if self._numWorkers > 1:
                    estimates = self._pool.map(ep_llda.run_vibound_EPLLDA,
                      [ (self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta, self._lambda,
                         self._omega, docs, words) for docs, words in zip(docchunks, wordchunks) ])
                else:
                    estimates = ep_llda.run_vibound_EPLLDA((self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta, self._lambda, self._omega, docchunks[0], wordchunks[0]))

            # "discriminative" variational
            elif (algname == 'discvar'):
                if self._numWorkers > 1:
                    res = self._pool.map(ep_llda.run_discvibound_EPLLDA,
                      [ (self._rngs[i], self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta,
                         self._lambda, self._omega, docs, words, self._Nsamp)
                         for i, docs, words in zip(range(numTasks),docchunks, wordchunks) ])
                    estimates = [ res[i][0] for i in range(numTasks) ]
                    self._rngs[0:numTasks] = [ res[i][1] for i in range(numTasks) ]
                else:
                    estimates, self._rngs[0] = ep_llda.run_discvibound_EPLLDA((self._rngs[0], self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta, self._lambda, self._omega, docchunks[0], wordchunks[0], self._Nsamp))

            # sample variational posterior
            elif (algname =='empirical'):
                if self._numWorkers > 1:
                    res = self._pool.map(ep_llda.run_estentropy_EPLLDA,
                      [ (self._rngs[i], self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta, self._lambda, self._omega, docs, words, self._Nsamp)
                        for i, docs, words in zip(range(numTasks), docchunks, wordchunks) ])
                    estimates = [ res[i][0] for i in range(numTasks) ]
                    self._rngs[0:numTasks] = [ res[i][1] for i in range(numTasks) ]
                else:
                    estimates, self._rngs[0] = ep_llda.run_estentropy_EPLLDA((self._rngs[0], self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta, self._lambda, self._omega, docchunks[0], wordchunks[0], self._Nsamp))

            # hybrid empirical / variational
            elif (algname == 'hybrid'):
                if self._numWorkers > 1:
                    res = self._pool.map(ep_llda.run_estMIhybrid_EPLLDA,
                      [ (self._rngs[i], self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta, self._lambda, self._omega, docs, words, self._Nsamp)
                        for i, docs, words in zip(range(numTasks), docchunks, wordchunks) ])
                    estimates = [ res[i][0] for i in range(numTasks) ]
                    self._rngs[0:numTasks] = [ res[i][1] for i in range(numTasks) ]
                else:
                    estimates, self._rngs[0] = ep_llda.run_estMIhybrid_EPLLDA((self._rngs[0], self._W, self._K, self._Nl, self._ppi, self._wordids, self._gamma, self._zeta, self._lambda, self._omega, docchunks[0], wordchunks[0], self._Nsamp))

            # find minimum conditional entropy
            MI = np.concatenate( estimates ).flatten()
            idx = np.argmax( MI )
            d_new, n_new = (doclist[idx], wordlist[idx])

        self._labels[d_new, n_new] = labels_all[d_new, n_new]
        return (d_new, n_new, self._labels[d_new, n_new], MI)


    def train(self):
        self.init_messages()

        # Prepare placeholder for old values of lambda and gamma
        old_lambda = np.empty_like(self._lambda)
        old_gamma = np.empty_like(self._gamma)

        # One iteration is going over all words in all documents
        for iteration in range(0,self._maxiters):
            # First make a copy of the old values to check for convergence later
            old_lambda[:] = self._lambda
            old_gamma[:] = self._gamma

            # Update omega and zeta. This method keeps lambda and gamma up to date with omega and zeta
            self.update_omega_zeta()

            # Compute new gamma
            delta = self.compute_diff(old_lambda, old_gamma)
            if not self._silent:
                print( "\tIter " + str(iteration+1) +"/" + str(self._maxiters) + \
                       ": delta = " + str(delta) + " (thr = " + str(self._threshold) + ")") #not pythonic

            # check convergence
            if delta < self._threshold:
                if not self._silent:
                    print("\tBatchEP converged in " + str(iteration+1) + \
                          " iterations, delta = " + str(delta) + " (thr = " + str(self._threshold) \
                          + ", avg_lambda = " + str(np.mean(self._lambda)) + ", avg_gamma = " \
                          + str(np.mean(self._gamma)) + ")")
                break
            if iteration >= self._maxiters - 1:
                print("\tBatchEP did not converge after " + str(iteration+1) + \
                      " iterations, delta = " + str(delta) + " (thr = " + str(self._threshold) \
                      + ", avg_lambda = " + str(np.mean(self._lambda)) + ", avg_gamma = " \
                      + str(np.mean(self._gamma)) + ")")
        return self._lambda

    def topic_entropy(self):
        """
        Compute entropy of topic posterior approx q(\phi)
        """
        Hhat = 0
        for k in range(self._K):
            Hhat += stats.dirichlet.entropy( self._lambda[k] )
        return Hhat

    def update_omega_zeta(self):
        """
        Iterates over all words in all documents to iteratively update
        omega and zeta, and keep lambda and gamma current.
        """
        for d in range(0, self._D):
            # print "\tAt doc " + str(d+1) + "/" + str(self._D) + ": N_d =  " + str(self._N[d])
            for n in range(0, self._N[d]):
                # print "\t\td = " + str(d+1) + "/" + str(self._D) + ", n = " + str(n+1) + "/" + str(self._N[d])

                # Compute rho_d = alpha + \sum_{n' \neq n} zeta_{dn'}
                # Equivalently, rho_d = gamma_d - zeta_dn
                rho_d = self._gamma[d,:] - self._zeta[d][n,:]  # array of length K

                # For 1 <= k <= K, compute tau_k = eta_k + \sum_{(d',n') \neq (d,n)} omega_{kd'n'}
                # Equivalently, tau_k = lambda_k - omega_kdn
                tau = np.empty_like(self._lambda)
                for k in range(0, self._K):
                    tau[k,:] = self._lambda[k,:] - self._omega[k][d][n,:]

                # If some components of rho_d or tau_k are non-positive, skip this step
                if (rho_d.min() <= 0 or tau.min() <= 0):
                    print("\t\tSkipped (d,n) = (" + str(d) + "," + str(n) \
                          + "): min(rho_d) = " + str(rho_d.min()) + ", min(tau) = " + str(tau.min()))
                    continue

                # Compute \sum_k \rho_dk and \sum_v \tau_kv, for 1 <= k <= K
                sum_rho_d = np.sum(rho_d)  # scalar
                sum_tau = np.sum(tau, 1)  # array of length K

                # The current word w_dn
                wdn = self._wordids[d][n]  # this is an integer between 0 and W-1 (size of the vocabulary)

                # Compute C_k, the (normalized) mixing proportions
                C = (rho_d * tau[:,wdn]) / sum_tau  # elementwise operation
                y_dn = self._labels[d,n]
                if (y_dn >= 0):
                    C *= self._ppi[:,y_dn]
                C = C / np.sum(C)  # nonnegative array of length K, summing to 1

                # Compute E[theta_dk] and E[theta_dk^2] with respect to the
                # approximate posterior distribution (which is a mixture of Dirichlet with
                # mixing proportion specified by C)
                E_theta_d = (rho_d + C) / (sum_rho_d + 1)
                E_theta_d2 = (rho_d * (rho_d + 1) + 2 * C * (rho_d + 1)) / ((sum_rho_d + 1) * (sum_rho_d + 2))

                # If we use Newton's method, we need to compute E[log theta_dk] as well
                if (self._useNewton):
                    E_log_theta_d = (1 - C) * psi(rho_d) + C * psi(rho_d + 1) - psi(sum_rho_d + 1)
                    self._gamma[d,:] = dirichlet_mle_newton(E_theta_d, E_theta_d2, E_log_theta_d)
                else:
                    self._gamma[d,:] = (sum(E_theta_d - E_theta_d2)) / (sum(E_theta_d2 - E_theta_d ** 2)) * E_theta_d

                # Then update zeta
                self._zeta[d][n,:] = self._gamma[d,:] - rho_d

                # For each k = 1, ..., K, compute E[beta_kv] and E[beta_kv^2] with respect to the
                # approximate posterior distribution, do moment-matching for beta_k, and update lambda_k and omega_k
                for k in range(0, self._K):
                    C_k = C[k]
                    tau_k = tau[k,:]
                    sum_tau_k = sum_tau[k]
                    tau_k_wdn = tau_k[wdn]

                    E_beta_k = tau_k / (sum_tau_k + 1) + (1 - C_k) * tau_k / (sum_tau_k * (sum_tau_k + 1))
                    E_beta_k[wdn] += C_k / (sum_tau_k + 1)

                    E_beta2_k = tau_k * (tau_k + 1) / ((sum_tau_k + 1) * (sum_tau_k + 2)) \
                                  + 2 * (1 - C_k) * tau_k * (tau_k + 1) / (sum_tau_k * (sum_tau_k + 1) * (sum_tau_k + 2))
                    E_beta2_k[wdn] += 2 * C_k * (tau_k_wdn + 1) / ((sum_tau_k + 1) * (sum_tau_k + 2))

                    # If we use Newton's method, we need to compute E[log beta_kv] as well
                    if (self._useNewton):
                        E_log_beta_k = psi(tau_k) - (1 - C_k) * psi(sum_tau_k) - C_k * psi(sum_tau_k + 1)
                        E_log_beta_k[wdn] += C_k * (psi(tau_k_wdn + 1) - psi(tau_k_wdn))
                        self._lambda[k,:] = dirichlet_mle_newton(E_beta_k, E_beta2_k, E_log_beta_k)
                    else:
                        self._lambda[k,:] = (sum(E_beta_k - E_beta2_k)) / (sum(E_beta2_k - E_beta_k ** 2)) * E_beta_k

                    # Then update omega
                    self._omega[k][d][n,:] = self._lambda[k,:] - tau[k,:]

    def  compute_diff(self, old_lambda, old_gamma):
        """
        Compute the average component-wise change
        in the values of current and old lambda and gamma.
        This method assumes self._lambda and self._gamma are up-to-date with
        the current values of self._omega and self._zeta.
        """
        gamma_diff = np.mean(abs(self._gamma - old_gamma))
        lambda_diff = np.mean(abs(self._lambda - old_lambda))
        return (0.5 * (gamma_diff + lambda_diff))

class LLDA_Gibbs:

    def __init__(self, pool, rngs, W, K, alpha, beta, ppi, Nl, labels, wordids, maxiters, \
                 threshold, numWorkers, useNewton, Nsamp, burn, \
                 phi_true, theta_true, z_true, silent=False):
        """
        Arguments:
        K: Number of topics
        W: Total number of words in the vocab
        alpha: Hyperparameter for prior on weight vectors theta
        beta: Hyperparameter for prior on topics beta
        Nl: Number of labels
        labels: Initial word labels
        wordids: list of lists with word indices
        maxiters: Max number of iterations to allow to converge
        threshold: Threshold for convergence
        numWorkers: Number of concurrent workers for parallel planning
        useNewton: Boolean, if true we use Newton's method for solving Dirichlet moment-matching, else use approximate MLE
        silent: (optional) Silence verbose output
        """

        self._str = "LLDA_EP_%d_%g_%r" % (maxiters,threshold,useNewton)
        self._K = K
        self._W = W
        self._alpha = alpha
        self._ppi = ppi
        self._Nl = Nl
        self._labels = labels
        self._maxiters = maxiters
        self._threshold = threshold
        self._numWorkers = numWorkers
        self._useNewton = useNewton
        self._silent = silent
        self._D = len(wordids)
        self._wordids = wordids
        self._zeta = list()
        self._omega = list()
        self._lambda = np.empty(0)
        self._gamma = np.empty(0)
        if np.isscalar(beta):
            self._beta = copy.deepcopy(beta) * np.ones((self._K, self._W))
        else:
            self._beta = copy.deepcopy(beta)

        # worker pool
        self._pool = pool
        self._rngs = rngs

        # true latents
        self._phi_true = phi_true
        self._theta_true = theta_true
        self._z_true = z_true

        # Get length of each document
        self._N = [len(ids) for ids in self._wordids]
        self._Nd = np.max( self._N )

        # sampler options
        self._Nsamp = Nsamp
        self._burn = burn

        # samples
        self._phi = np.zeros((self._K,self._W,self._Nsamp))
        self._theta = np.zeros((self._D,self._K,self._Nsamp))
        self._z = np.zeros((self._Nd,self._D,self._Nsamp), dtype='int')
        self._count_z = np.zeros((self._W,self._K,self._Nsamp), dtype='int')

    def train(self):
        """
        Run independent parallel Gibbs chains
        """

        debug = False

        # fit model in parallel
        if self._numWorkers > 1:
            this_Nsamp = int(math.ceil(self._Nsamp / float(self._numWorkers)))
            res = self._pool.map(gibbs_llda.gibbs,
                [ (self._rngs[i], self._W, self._K, self._D, self._Nd, self._Nl, self._labels,
                   self._wordids, self._alpha, self._beta, self._ppi, this_Nsamp, self._burn)
                for i in range(self._numWorkers) ] )

            # unpack results
            phi_map_all = np.concatenate([ res[i][0][...,np.newaxis] for i in range(self._numWorkers) ], axis=2)
            logp_map_all = [ res[i][1] for i in range(self._numWorkers) ]
            self._phi = np.concatenate([ res[i][2] for i in range(self._numWorkers) ], axis=2)
            self._theta = np.concatenate([ res[i][3] for i in range(self._numWorkers) ], axis=2)
            self._z = np.concatenate([ res[i][4] for i in range(self._numWorkers) ], axis=2)
            self._count_z = np.concatenate([ res[i][5] for i in range(self._numWorkers) ], axis=2)
            logp_trace = np.concatenate([ res[i][6] for i in range(self._numWorkers) ], axis=0)
            self._rngs = [ res[i][7] for i in range(self._numWorkers) ]

            # throw away extra samples
            self._phi = self._phi[:,:,0:self._Nsamp]
            self._theta = self._theta[:,:,0:self._Nsamp]
            self._z = self._z[:,:,0:self._Nsamp]
            self._count_z = self._count_z[:,:,0:self._Nsamp]
            logp_trace = logp_trace[0:self._Nsamp,:]

            # compute MAP of MAPs
            idx_max = np.argmax( logp_map_all )
            phi_map = phi_map_all[...,idx_max]
            logp_map = logp_map_all[idx_max]

        # single-threaded fitting
        else:
            phi_map, logp_map, self._phi, self._theta, self._z, self._count_z, \
                logp_trace, self._rngs[0] = gibbs_llda.gibbs( (self._rngs[0], self._W,
                self._K, self._D, self._Nd, self._Nl, self._labels, self._wordids,
                self._alpha, self._beta, self._ppi, self._Nsamp, self._burn))

        # DEBUG: Show traceplots
        if debug:
            plt.figure()
            for isamp in range(self._Nsamp):
                plt.plot(range(self._burn), logp_trace[isamp,:], '-b')
            plt.xlabel('Gibbs Iteration')
            plt.ylabel('Log-Joint Probability')

            # # compute true probability
            # logp_true = gibbs_llda.eval_logjoint(self._alpha, self._beta, self._ppi, self._labels,
            #                                 self._wordids, self._phi_true, self._theta_true)
            # plt.plot(range(self._burn), logp_true * np.ones(self._burn), '-r')

            # DEBUG: Show TV error vs. probability
            err = np.zeros(self._Nsamp)
            for isamp in range(self._Nsamp):
                err[isamp] = utils.tv_error( self._phi_true,  self._phi[:,:,isamp] )
            err_map = utils.tv_error(self._phi_true, phi_map)
            plt.figure()
            plt.scatter(err, logp_trace[:,-1])
            plt.scatter(err_map, logp_map, label='MAP')
            plt.xlabel('Topic Error (TV)')
            plt.ylabel('Log-Probability')

            # show MAP topics
            utils.show_bars_topics(phi_map, self._W, self._K)
            plt.show()

        # # DEBUG: Show all samples *before* relable
        # if debug:
        #     for isamp in range(self._Nsamp):
        #         fig, axs = utils.show_bars_topics(self._phi[...,isamp], self._W, self._K)
        #         fig.suptitle('Sample %d (Before)' % isamp)

        # return posterior mean estimate
        self.relabel_samples(self._phi_true)
        phi_avg = 1 / self._Nsamp * np.sum(self._phi, axis=2)

        # # DEBUG: Show all samples *after* relable
        # if debug:
        #     fig, axs = utils.show_bars_topics(phi_avg, self._W, self._K)
        #     fig.suptitle('Average')
        #     fig, axs = utils.show_bars_topics(phi_map, self._W, self._K)
        #     fig.suptitle('MAP')
        #     for isamp in range(self._Nsamp):
        #         fig, axs = utils.show_bars_topics(self._phi[...,isamp], self._W, self._K)
        #         fig.suptitle('Sample %d (After)' % isamp)
        #     plt.show()

        return phi_avg

    def relabel_samples(self, phi_target):
        """
        Relabel topics on all samples to be consistent with target.
        """

        sizeOfChunks = int(math.ceil(self._Nsamp / float(self._numWorkers)))
        phi_chunks = utils.chunk_array(self._phi, sizeOfChunks)
        theta_chunks = utils.chunk_array(self._theta, sizeOfChunks)
        count_z_chunks = utils.chunk_array(self._count_z, sizeOfChunks)
        z_chunks = utils.chunk_array(self._z, sizeOfChunks)
        numTasks = len(phi_chunks)

        if self._numWorkers > 1:
            res = self._pool.map(gibbs_llda.relabel_samples,
                [ (phi_target, phi_chunks[i], theta_chunks[i], count_z_chunks[i], z_chunks[i])
                    for i in range(numTasks) ] )

            self._phi = np.concatenate([ res[i][0] for i in range(numTasks) ], axis=2)
            self._theta = np.concatenate([ res[i][1] for i in range(numTasks) ], axis=2)
            self._count_z = np.concatenate([ res[i][2] for i in range(numTasks) ], axis=2)
            self._z = np.concatenate([ res[i][3] for i in range(numTasks) ], axis=2)

        else:
            self._phi, self._theta, self._count_z, self._z = gibbs_llda.relabel_samples(
                ( phi_target, phi_chunks[0], theta_chunks[0], count_z_chunks[0], z_chunks[0]) )


    def label_element(self, labels_all, algname):
        """
        Empirical estimate of conditional entropy H(\psi | Y)
        """

        # find unlabeled words and "chunk" them
        doclist, wordlist = np.where( self._labels == -1 )
        sizeOfChunks = int(math.ceil(len(doclist) / float(self._numWorkers)))
        docchunks = utils.chunk(doclist, sizeOfChunks)
        wordchunks = utils.chunk(wordlist, sizeOfChunks)
        numTasks = len(docchunks)

        # random planning
        if (algname == 'random'):
            N = len( doclist )
            idx = np.random.randint( N )
            d_new, n_new = (doclist[idx], wordlist[idx])
            Hhat = []

        # build argument lists for each worker
        else:

            # parallel planning
            if self._numWorkers > 1:
                args = [(self._rngs[i], self._phi, self._z, self._count_z, self._ppi, self._K,
                    self._Nl, self._Nd, self._beta, docs, words)
                    for i, docs, words in zip(range(numTasks), docchunks, wordchunks)]
                res = self._pool.map(gibbs_llda.run_estEntropyGibbs, args)

                # unpack retvals
                Hhat = [ res[i][0] for i in range(numTasks) ]
                self._rngs[0:numTasks] = [ res[i][1] for i in range(numTasks) ]

            # sequential planning
            else:
                Hhat, self._rngs[0] = gibbs_llda.run_estEntropyGibbs((self._rngs[0], self._phi, self._z,
                    self._count_z, self._ppi, self._K, self._Nl, self._Nd, self._beta,
                    docchunks[0], wordchunks[0]))

            # estimate etropy
            Hhat = np.concatenate(Hhat).flatten()

            # rank conditional entropy
            idx = np.argmin(Hhat)
            d_new, n_new = (doclist[idx], wordlist[idx])

        # return stuff
        self._labels[d_new, n_new] = labels_all[d_new, n_new]
        return (d_new, n_new, self._labels[d_new, n_new], Hhat)

    def topic_entropy(self):
        """
        Estimate entropy of topic distribution on current samples.
        """

        # split samples
        Nsamp = int(np.ceil(self._Nsamp / 2)) # self._Nsamp #
        count_z_samp_marg = self._count_z[:, :, Nsamp:]  # W x K x Nsamp # self._count_z #
        phi_samp = self._phi[:, :, 0:Nsamp] # K, W, Nsamp # self._phi #

        # compute Dirichlet conditional p(\phi | z, words)
        log_pcond_phi = np.zeros((self._K,Nsamp, Nsamp))  # K x phi-samples x z-samples
        beta_sim = self._beta.T[:, :, np.newaxis] + count_z_samp_marg  # W x K x Nsamp
        for k in range(self._K):
            log_pcond_phi[k] = \
                utils.dirichlet_logpdf(phi_samp[k, ...], beta_sim[:, k, :])

        # estimate marginal p(\phi | words)
        logZ = np.max( log_pcond_phi, axis=2 ) # K x phi-samples
        pcond_phi = np.exp( log_pcond_phi - logZ[...,np.newaxis] ) # K x phi-samples x z-samples
        logp_phi = logZ + np.log( 1/Nsamp * np.sum(pcond_phi, axis=2) ) # K x phi-samples

        # estimate entropy
        Hhat = 1/Nsamp * np.sum( - logp_phi, axis=1 )
        return Hhat
