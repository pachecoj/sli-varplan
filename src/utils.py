# utils.py: Utility functions

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
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details. You should have received a copy of the GNU General
# Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.

import sys, re, time, string
import copy
import numpy as n
from scipy.special import psi, polygamma, gammaln
from scipy.stats import stats
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt

# This is used in BatchVB's do_e_step
meanchangethresh = 0.00001

def unzipDocs(docs):
    wordids = [ids for (ids,cts) in docs]
    wordcts = [cts for (ids,cts) in docs]
    return (wordids,wordcts)


def dirichlet_expectation(alpha):
    """
    For a vector theta ~ Dir(alpha), computes E[log theta] and E[theta] given alpha.
    """
    if (len(alpha.shape) == 1):
        return( (psi(alpha) - psi(n.sum(alpha))),
                alpha / n.sum(alpha) )
    return ( (psi(alpha) - psi(n.sum(alpha, 1))[:, n.newaxis]),
             ( alpha / n.sum(alpha, 1)[:, n.newaxis] ) )

def strToBool(s):
    """
    Convert a string s into a Boolean value.
    """
    return (s.lower() in ['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh', 'indubitably'])


def unzipDocsShuffle(docs, labels):
    """
    Extracts the (shuffled) list of wordids (i.e. the integer ids from the
    vocab list) for each doc in docs. The output format is a list of lists,
    one for each document. Note that the length may be different for each
    document.
    """
    D = len(docs)  # number of documents
    wordids = list()  # placeholder for output
    labels_new = n.zeros( labels.shape )

    for d in range(0, D):
        (ids, cts) = docs[d]
        Nd = sum(cts)  # number of words in doc d

        # First populate wordids_d with the ordered ids,
        # e.g. if ids = [3,8,20] and cts = [1,3,2], then
        # wordids_d = [3,8,8,8,20,20]
        wordids_d = list()
        for i in range(0, len(ids)):
            wordids_d.extend([ids[i]] * cts[i])

        # Shuffle wordids_d and insert to wordids
        idx = n.random.permutation(range(len(wordids_d)))
        wordids_d = [ wordids_d[i] for i in idx ]
        labels_new[d,:] = labels[d,idx]
        wordids.append(wordids_d)
    return (wordids, labels_new)



def dirichlet_mle_newton(e_p, e_p2, e_logp, maxiters = 20, thr = 1e-4, silent = False):
    """
    Finds the MLE for the K-dimensional Dirichlet distribution from observed data,
    i.e. the solution alpha_1, ..., alpha_K > 0 to the moment-matching equations
        psi(alpha_k) - psi(sum(alpha)) = E[log p_k]
    where the expectation on the right hand side is with respect to the empirical
    distribution.

    Input: e_p, a vector of length K containing the empirical expectations E[p_k], i.e. e_p.ndim == 1 and len(e_p) == K
           e_p2, the empirical expectations E[p_k^2], the same format as e_p
           e_logp, the empirical expectations E[log p_k], the same format as e_p
           maxiters, the maximum number of Newton-Raphson iterations
           thr, the threshold for convergence
    Output: alpha, a vector of length K containing the parameters alpha_1, ..., alpha_K

    This method uses the first and second empirical moments e_p and e_p2 to initialize
    the alpha values (by approximately matching the first and second moments), and then
    uses Newton-Raphson method to refine the estimates.

    This method is based on the first section of Minka's paper:
    http://research.microsoft.com/en-us/um/people/minka/papers/dirichlet/minka-dirichlet.pdf
    """

    # For initialization: First compute the approximate sum(alpha)
    alpha0 = (sum(e_p - e_p2)) / (sum(e_p2 - e_p ** 2))

    # Then compute the initial alpha
    alpha = alpha0 * e_p

    # Do Newton-Raphson iterations
    for iteration in range(0, maxiters):
        sum_alpha = sum(alpha)
        g = psi(alpha) - psi(sum_alpha) - e_logp
        z = polygamma(1, sum_alpha)  # polygamma(1,z) is the trigamma function psi_1(z)
        q = polygamma(1, alpha)
        b = sum(g / q) / (1 / z - sum(1 / q))
        alpha_new = alpha - (g + b) / q

        # this is a hack, but if some of alpha_new's components are negative, make them positive
        # alpha_new[alpha_new < 0] = alpha / 5  # / 5 is arbitrary, as long as the end result is positive

        # Update alpha and check for convergence
        delta = max(abs(alpha - alpha_new))
        alpha = alpha_new
        if delta < thr:
            # cur_gap = psi(alpha) - psi(sum(alpha)) - e_logp
            # if not silent:
            #     print "Dirichlet-MLE-Newton converged in " + str(iteration) + " iterations, gap = " + str(cur_gap)
            break
        if iteration >= maxiters - 1:
            cur_gap = psi(alpha) - psi(sum(alpha)) - e_logp
            if not silent:
                print("Dirichlet-MLE-Newton did not converge after " + \
                      str(iteration) + " iterations, gap = " + str(cur_gap))
    return alpha


def sample_labels(labels_all, Nsel, Nlabels):
    """
    Samples N labels from 2D Numpy array labels_all.  Returns a 2D
    Numpy array of selected labels and NaN's for unselected labels.
    """
    N, M = labels_all.shape
    IJ = n.random.randint(Nlabels, size=(N,M))
    idx = n.random.permutation(N*M)
    idx = idx[0:Nsel]
    I, J = n.unravel_index(idx, (N,M), order='F')
    labels = int(-1) * n.ones((N,M), dtype=int)
    labels[ I, J ] = labels_all[ I, J ]
    return labels, list(I), list(J)


def show_bars_topics(lam, W, K):
    """
    Display "bars" topics in subplot figure.
    """
    fig_top, axs_top = plt.subplots(2,int(K/2))
    axs_top = axs_top.reshape(K)
    for k in range(K):
        this_lam = n.reshape( lam[k,:], ( int( n.sqrt(W) ), int( n.sqrt(W) ) ) ) / n.sum( lam[k,:] )
        axs_top[k].imshow( this_lam, interpolation="none", cmap="gray" )
        axs_top[k].set_title('Topic %d' % k)
        axs_top[k].set_xticks([])
        axs_top[k].set_yticks([])
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    return (fig_top, axs_top)

def tv_error(phi_true, phi_est):
    """
    Compute total variation error and account for label switching
    """
    tmp1, tmp2, tv = align_topics(phi_true, phi_est)
    return tv

def align_topics(phi_true, phi):
    """
    Solve for optimal ordering of topics by minimizing total variation
    """
    K, W = phi_true.shape

    # compute costs
    cost = n.zeros((K,K))
    for k1 in range(K):
        top_true = phi_true[k1,:]
        for k2 in range(K):
            top = phi[k2,:]
            cost[k1,k2] = 1/2 * n.sum( n.abs( top - top_true ) )

    # do matching
    row, col = linear_sum_assignment(cost)
    phi_new = phi[col,:]
    tv_err = cost[row,col].sum() / K
    return (phi_new, col, tv_err)


def topic_entropy_vec(lam, K):
    H = n.zeros((K,))
    for k in range(K):
        H[k] = stats.dirichlet.entropy(lam)
    return H

def dirichlet_logpdf(X, a):
    """
    Compute log-pdf of a Dirichlet distribution
    Input:
      X: d x n data matrix, each column sums to one (sum(X,1)==ones(1,n) && X>=0)
      a: d x k parameter of Dirichlet
    Output:
      y: k x n probability density in logrithm scale y=log p(x)
    Adapted from Mo Chen (sth4nth@gmail.com)
    """
    if a.ndim == 1:
        a = a[:,n.newaxis]
    c = gammaln(n.sum(a,axis=0))-n.sum(gammaln(a),axis=0) # k-vec or scalar
    g = n.dot((a.T-1), n.log(X)) # k x n
    y = g + c[:,n.newaxis] # k x n
    return y
