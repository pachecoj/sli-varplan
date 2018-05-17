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
corpus = 'TINYBARS'
batchsize = 1024
algorithmname = 'llda_ep'  # 'llda_ep' or 'llda_gibbs'
sel_algorithmname = 'variational'  # 'random' or 'variational' or 'empirical' or 'hybrid'
code = 'REALENT'
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
    err = np.zeros((train_steps,))
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
        print('TV Error: %0.2f' % err[t])
        sio.savemat(this_resfile, {'Hcond':Hcond, 'Hreal':Hreal, 'err':err, 'lambda':lam, 'W':W,
            'K':K, 'max_iters':max_iters, 'threshold':threshold, 'useNewton':useNewton, 'sel_d':sel_d, 'sel_n':sel_n, 'corpus':corpus})

    # done
    print('Saved %s' % this_resfile)

# #     # show topics
# # utils.show_bars_topics(alg_lam, W, K)
# # plt.figure()
# # plt.plot( range(train_steps), err )
# # plt.xlabel('Training Epoch')
# # plt.title('Topic TV Error')
# # plt.show()


#### SCRAP ####################

# if (sel_algorithmname == 'random'):
            #     d_new, n_new, l_new = alg.label_element_random(labels_all)
            #     sel_d.append( d_new )
            #     sel_n.append( n_new )
            #
            # elif (sel_algorithmname == 'variational'):
            #
            #     labels_old = alg._labels.copy()  # for debugging
            #     d_new, n_new, l_new, this_H = alg.label_element_vi(labels_all, pool)
            #     Hcond.append( this_H )
            #     sel_d.append( d_new )
            #     sel_n.append( n_new )
            #
            #     ## DEBUG: Calculate hybrid estimate
            #     labels_new = alg._labels.copy()
            #     MI_hyb = []
            #     Nmc = 10
            #     for i in range(Nmc):
            #         print('Sampling entropy %d' % i)
            #         alg._labels = labels_old.copy()
            #         tmp1, tmp2, tmp3, this_MI = alg.label_element_hybrid(Nsamp, labels_all, pool)
            #         MI_hyb.append(this_MI)
            #     MI_hyb = np.array(MI_hyb)
            #     MI_mean_hyb = np.mean(MI_hyb, axis=0)
            #     MI_std_hyb = np.std(MI_hyb, axis=0)
            #
            #     ## DEBUG: Calculate empirical estimate
            #     MI_emp = []
            #     for i in range(Nmc):
            #         print('Sampling entropy %d' % i)
            #         alg._labels = labels_old.copy()
            #         tmp1, tmp2, tmp3, this_MI = alg.label_element_emp(Nsamp, labels_all, pool)
            #         MI_emp.append(this_MI)
            #     MI_emp = np.array(MI_emp)
            #     MI_mean_emp = np.mean(MI_emp, axis=0)
            #     MI_std_emp = np.std(MI_emp, axis=0)
            #     alg._labels = labels_new.copy()
            #
            #     ## DEBUG: Plot stuff
            #     plt.figure()
            #     plt.plot(range(len(this_MI)), MI_mean_emp, '-r', label='Empirical')
            #     plt.plot(range(len(this_MI)), MI_mean_emp - MI_std_emp, '--r', label='STDEV')
            #     plt.plot(range(len(this_MI)), MI_mean_emp + MI_std_emp, '--r')
            #     plt.plot(range(len(this_MI)), MI_mean_hyb, '-b', label='Hybrid')
            #     plt.plot(range(len(this_MI)), MI_mean_hyb - MI_std_hyb, '--b', label='STDEV')
            #     plt.plot(range(len(this_MI)), MI_mean_hyb + MI_std_hyb, '--b')
            #     plt.plot(range(len(this_MI)), this_H, '-k', label='Full Variational')
            #     plt.xlabel('Word Index')
            #     plt.ylabel('Mutual Information')
            #     plt.title('Epoch %d, Samples %d' % (t,Nmc))
            #     plt.legend()
            #     plt.show()
            #
            # elif (sel_algorithmname == 'empirical'):
            #     d_new, n_new, l_new, this_H = alg.label_element_emp(Nsamp, labels_all, pool)
            #     Hcond.append( this_H )
            #     sel_d.append( d_new )
            #     sel_n.append( n_new )
            #
            # elif (sel_algorithmname == 'hybrid'):
            #     labels_old = alg._labels.copy() # for debugging
            #     d_new, n_new, l_new, this_H = alg.label_element_hybrid(Nsamp, labels_all, pool)
            #     Hcond.append( this_H )
            #     sel_d.append( d_new )
            #     sel_n.append( n_new )
            #
            #     # ## DEBUG: Calculate hybrid estimate
            #     # labels_new = alg._labels.copy()
            #     # MI_emp = []
            #     # Nmc = 10
            #     # for i in range(Nmc):
            #     #     print('Sampling entropy %d' % i)
            #     #     alg._labels = labels_old.copy()
            #     #     tmp1, tmp2, tmp3, this_MI = alg.label_element_emp(Nsamp, labels_all, pool)
            #     #     MI_emp.append(this_MI)
            #     # MI_emp = np.array(MI_emp)
            #     # MI_mean = np.mean(MI_emp, axis=0)
            #     # MI_std = np.std(MI_emp, axis=0)
            #     # alg._labels = labels_new.copy()
            #     #
            #     # ## DEBUG: Plot stuff
            #     # plt.figure()
            #     # plt.plot(range(len(this_MI)), MI_mean, '-b', label='Empirical')
            #     # plt.plot(range(len(this_MI)), MI_mean - MI_std, '--b', label='STDEV')
            #     # plt.plot(range(len(this_MI)), MI_mean + MI_std, '--b')
            #     # plt.plot(range(len(this_MI)), this_H, label='Hybrid')
            #     # plt.xlabel('Word Index')
            #     # plt.title('MI (Epoch %d, Samples %d)' % (t,Nmc))
            #     # plt.legend()
            #     # plt.show()
            #
            # else:
            #     # raise Exception("Selection algorithm %s not recognized." % sel_algorithmname)


            # # DEBUG: Show selection
            # wdn = training_wordids[d_new][n_new]
            # fig, axs = utils.show_bars_topics(lam[:,:,t-1], W, K)
            # row, col = np.unravel_index(wdn, (5,5))
            # axs[l_new].plot(col, row, '*r')
            # fig.suptitle('Epoch %d Selection (%s)' % (t,sel_algorithmname))
            # plt.savefig('/home/pachecoj/Dropbox/Research/journal/042218/tmp/sel%d.png' % t)
            # plt.close()

            # # DEBUG: Show entropies
            # import scipy.stats as stats
            # labels_orig = alg._labels.copy()
            # d_all, n_all = np.where(labels_orig == -1)
            # H_avg = np.ones((K,W))
            # C = np.ones((K,W))
            # for i in range(len(d_all)):
            #     w = training_wordids[d_all[i]][n_all[i]]
            #     k = int( labels_all[d_all[i], n_all[i]] )
            #     H_avg[k,w] += this_H[i]
            #     C[k,w] += 1
            # H_avg = H_avg / C
            # fig, axs = utils.show_bars_topics(H_avg, W, K)
            # wdn = training_wordids[d_new][n_new]
            # row, col = np.unravel_index(wdn, (5,5))
            # axs[l_new].plot(col, row, '*r')
            # fig.suptitle('Entropy Epoch %d (%s)' % (t, sel_algorithmname))
            # plt.savefig('/afs/csail.mit.edu/u/p/pachecoj/Dropbox/Research/journal/042218/sel_%s/H%d.png' % (sel_algorithmname,t))
            # plt.close()

            # # DEBUG: Show entropy of theta_dk
            # import scipy.stats as stats
            # import time
            # labels_orig = alg._labels.copy()
            # d_all, n_all = np.where(labels_orig == -1)
            # d_all = np.unique( d_all )
            # tmp, theta_est = utils.dirichlet_expectation(alg._gamma)
            # for d in d_all:
            #     plt.bar(range(K), theta_est[d,:])
            #     plt.title('Doc %d' % d)
            #     plt.xlabel('Topic')
            #     plt.ylabel('Doc/Topic Entropy')
            #     plt.show()
            #     time.sleep(1.)
