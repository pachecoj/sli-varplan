import pickle
import batch, utils
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
from llda_model import *
from multiprocessing import Pool, cpu_count
import archived_dataset

np.random.seed(100000001)

# set options
datadir = '../data/'
corpus = 'SPARSEBARS'
algorithmname = 'llda_ep'  # 'llda_ep' or 'llda_gibbs'
sel_algorithmname = 'random'  # 'discvar' or 'random' or 'variational' or 'empirical' or 'hybrid'
code = 'DEBUG_DISCVAR'
numinit = 50
max_iters = 100
threshold = 0.1
train_steps = 100
useNewton = False
Nsamp = 1000
burn = 100
numWorkers = cpu_count()
Nruns = 10
pool = Pool(numWorkers)

# load vocab / training / test files
W = len(open(datadir + corpus + "_vocab.dat", 'r').readlines())
training_docs = archived_dataset.Corpus(datadir + corpus + "_train.dat")
# validation_docs = archived_dataset.loadDocs(datadir + corpus + "_test.dat")

# load model
f = open(datadir + corpus + '_model.dat', 'rb')
model = pickle.load(f)
f.close()
eta = model.beta     # prior for topics
K = model.K          # number of topics
alpha = model.alpha  # topic proportions prior
ppi = model.ppi      # label conditional
Nl = model.Nl        # number of labels
D = training_docs._D # number of documents

# results filename
resfile = '%scorpus_%s_alg_%s_sel_%s_numinit_%d_trainsteps_%d_Nsamp_%d_%s.mat' % \
    (datadir, corpus, algorithmname, sel_algorithmname, numinit, train_steps, Nsamp, code)

# run a bunch of times
for run_id in range(Nruns):
    this_resfile = '%s_runid_%d.mat' % (resfile[0:-4], run_id)

    # load all labels and shuffle
    labels_all = archived_dataset.loadLabels(\
                    datadir + corpus + "_labels.dat")
    training_wordids, labels_all = utils.unzipDocsShuffle(training_docs._data, labels_all)

    # pre-generate random states for Pool tasks
    maxTasks = 10 * numWorkers # shouldn't be more than this
    rngs = [np.random.RandomState(run_id*maxTasks + i) for i in range(maxTasks)]

    # sample labels
    if numinit:
        labels, sel_d, sel_n = utils.sample_labels(labels_all, numinit, model.Nl)

    # # DEBUG: Init all but topics 2 & 7
    # I, J = np.where(( labels_all != 2 ) & (labels_all != 7))
    # labels = int(-1) + np.zeros_like(labels_all, dtype=int)
    # labels[I,J] = labels_all[I,J]
    # numinit = len(I)

    # # DEBUG: Init all but first word
    # labels = labels_all.copy()
    # labels[0,0] = -1
    # numinit = len( labels.flatten() ) - 1

    # init algorithm
    if (algorithmname == "llda_ep"):     # EP for Labeled LDA
        alg = batch.LLDA_EP(
                pool, rngs, W, K, alpha, eta, ppi, model.Nl, labels, training_wordids, max_iters,
                threshold, numWorkers, useNewton, Nsamp, silent=False)
    elif (algorithmname == "llda_gibbs"):     # Gibbs for Labeled LDA
        alg = batch.LLDA_Gibbs(
                pool, rngs, W, K, alpha, eta, ppi, model.Nl, labels, training_wordids, max_iters,
                threshold, numWorkers, useNewton, Nsamp, burn,
                model.phi_true, model.theta_true, model.z_true, silent=False)
    else:
        raise Exception("Algorithm %s not recognized." % algorithmname)

    # fit training data
    err = np.nan + np.zeros((train_steps,))
    Hcond = []
    Hreal = []
    lam = np.zeros((K,W,train_steps))
    for t in range(train_steps):
        print('Training Epoch %d' % t)

        # planning
        if t>0:
            d_new, n_new, l_new, this_H = alg.label_element(labels_all, sel_algorithmname)
            Hcond.append( this_H )
            sel_d.append( d_new )
            sel_n.append( n_new )

        # inference
        lam[:,:,t] = alg.train()
        tmp, phi_est = utils.dirichlet_expectation(lam[:,:,t])
        err[t] = utils.tv_error( model.phi_true,  phi_est )

        # compute posterior entropy
        Hreal.append( alg.topic_entropy() )

        # save results
        print('TV Error: %0.4f' % err[t])
        print('Realized Entropy: %0.2f' % np.sum(Hreal[-1]))
        sio.savemat(this_resfile, {'Hcond':Hcond, 'Hreal':Hreal, 'err':err, 'lambda':lam, 'W':W,
            'K':K, 'max_iters':max_iters, 'threshold':threshold, 'useNewton':useNewton, 'sel_d':sel_d, 'sel_n':sel_n, 'corpus':corpus})

    # done
    print('Saved %s' % this_resfile)
