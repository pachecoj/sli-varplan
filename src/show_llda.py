import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import utils

datadir = '../data/'

files = ('corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_10_trainsteps_10_Nsamp_100_SHORT_runid_0.mat',
         'corpus_SPARSEBARS_alg_llda_ep_sel_discvar_numinit_10_trainsteps_10_Nsamp_100_SHORT_runid_0.mat')
names = ('Gibbs','Variational',)
colors = ('-b','-k',)
show_iters = np.arange(0,10)

# plot results
fig_tv = plt.figure()
for fname, algname, idx in zip(files, names, range(len(files))):
    print(datadir + fname)
    res = sio.loadmat( datadir + fname )
    err = res['err'].flatten()
    plt.figure( fig_tv.number )
    plt.plot( range(len(err)), err, colors[idx], label=algname)

    # plot final topics
    this_lambda = res['lambda']    
    fig_top, ax_top = utils.show_bars_topics(this_lambda[:,:,-1], int(res['W'][0]), int(res['K'][0]))
    fig_top.suptitle('%s Topics' % algname)
    
    # # plot topics
    # for i in show_iters:        
    #     fig_top, ax_top = utils.show_bars_topics(this_lambda[:,:,i], int(res['W'][0]), int(res['K'][0]))
    #     fig_top.suptitle('%s Topics (Iter %d)' % (algname, i))

    #     # # show only 1st topic for sparse bars
    #     # fig_top = utils.show_single_bars_topic(this_lambda[...,i], int(res['W'][0]), 0)
    #     # fig_top.suptitle('%s Topic (Iter %d)' % (algname, i))
    #     # # fig_top.savefig('/home/pachecoj/Dropbox/Research/journal/091618/%s_iter%d.pdf' % (algname,i))

plt.figure( fig_tv.number )
plt.legend()
plt.xlabel('Training Epoch')
plt.ylabel('Topic Error (TV)')
plt.show()
