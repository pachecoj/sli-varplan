import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import pickle
import utils, archived_dataset

datadir = '../data/'
fname = 'corpus_TINYBARS_alg_llda_batch_ep_sel_var_numinit_10_trainsteps_100_res.mat'

# load results
res = sio.loadmat( datadir + fname )

# load model / data
corpus = 'TINYBARS'
# f = open(datadir + corpus + '_model.dat', 'rb')
# model = pickle.load(f)
# f.close()
training_docs = archived_dataset.Corpus(datadir + corpus + "_train.dat")

# unpack stuff
# Note:  Assumes all docs have same number of words!
# TODO: Load these from data
err = res['err'].flatten()
D = len( training_docs._data )
Nd = 25
numinit = 10
train_steps = len( res['Hcond'][0,:] )

# get selected indices
sel_d = res['sel_d'][0]
sel_n = res['sel_n'][0]
sel_idx = np.ravel_multi_index((sel_d, sel_n), (D,Nd), order='F')

# plt.figure( fig_tv.number )
# plt.plot( range(len(err)), err)
#
# # plot topics
# this_lambda = res['lambda']
# fig_top, ax_top = utils.show_bars_topics(this_lambda[:,:,-1], int(res['W'][0]), int(res['K'][0]))
# fig_top.suptitle('%s Topics Variational')
# print( res['sel'] )

# plot entropy
plt.figure()
Hcond = res['Hcond'][0,:]
for t in range(train_steps):
    print('sel: (%d, %d)' % (sel_d[t], sel_n[t]))

    # get indices
    idx_used = sel_idx[0:(numinit + t)]
    idx_avail = np.arange(D * Nd)
    # if t == 28:
    #     import ipdb
    #     ipdb.set_trace()
    idx_avail = np.delete(idx_avail, idx_used)

    # set conditional entropy values
    this_Hcond = np.empty((D*Nd,))
    this_Hcond[:] = np.NaN
    this_Hcond[idx_avail] = Hcond[t].flatten()

    # plot
    plt.plot(range(D*Nd), this_Hcond)

plt.xlabel('Word Index')
plt.ylabel('Conditional Entropy H(\phi | Y)')

# plot histogram of selections
hist = np.histogram2d(res['sel_d'][0], res['sel_n'][0], bins=(np.arange(D), np.arange(Nd)))
plt.figure()
plt.imshow( hist[0], interpolation='none', cmap='gray' )
plt.xlabel('Word Index')
plt.ylabel('Document Index')
plt.title('Variational Selections')

# show plot
plt.show()
