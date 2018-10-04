import matplotlib.pyplot as plt
import matplotlib.colors as pycolor
import scipy.io as sio
import numpy as np
import utils

datadir = '../data/'

files = ('corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_150_trainsteps_100_Nsamp_100_COUNT_SEL_runid_0.mat',
         'corpus_SPARSEBARS_alg_llda_ep_sel_discvar_numinit_150_trainsteps_100_Nsamp_100_COUNT_SEL_runid_0.mat',)
names = ('Gibbs','Variational',)
colors = ('-b','-k',)
show_iters = np.arange(240,251)

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

    # get selected words
    corpus = np.array( res['training_wordids'] )
    d = res['sel_d']
    n = res['sel_n']
    W = int( res['W'][0] )

    # get the scale factor for selection
    start_iter = show_iters[0]
    stop_iter = show_iters[-1]
    sel_wordids = corpus[d[0,start_iter:stop_iter],n[0,start_iter:stop_iter]]
    sel_vec, trash = np.histogram(sel_wordids, bins=np.arange(0,W+1))
    scale = np.max( sel_vec )
        
    # plot topics
    for i in show_iters:

        # histogram selected words
        sel_wordids = corpus[d[0,start_iter:i],n[0,start_iter:i]]
        sel_vec, trash = np.histogram(sel_wordids, bins=np.arange(0,W+1))
        f_selvec = sel_vec / float( scale )
        fig_top = plt.figure()
        this_lam = np.reshape( f_selvec, ( int( np.sqrt(W) ), int( np.sqrt(W) ) ) )
        plt.imshow( this_lam, interpolation="none", cmap="gray", vmin=0., vmax=1. )            
        # fig_top = utils.show_single_bars_topic(fsel_vec[np.newaxis,:], W, 0)
        fig_top.suptitle('%s Selections (Iter %d)' % (algname, i))
        fig_top.savefig('/home/pachecoj/Dropbox/Research/journal/093018/%s_selection_iter%d.pdf' % (algname,i))
        plt.close( fig_top )
        
        # fig_top, ax_top = utils.show_bars_topics(this_lambda[:,:,i], int(res['W'][0]), int(res['K'][0]))
        # fig_top.suptitle('%s Topics (Iter %d)' % (algname, i))

        # # show only 1st topic for sparse bars
        # fig_top = utils.show_single_bars_topic(this_lambda[...,i], int(res['W'][0]), 0)
        # fig_top.suptitle('%s Topic (Iter %d)' % (algname, i))
        # # fig_top.savefig('/home/pachecoj/Dropbox/Research/journal/091618/%s_iter%d.pdf' % (algname,i))

plt.figure( fig_tv.number )
plt.legend()
plt.xlabel('Training Epoch')
plt.ylabel('Topic Error (TV)')
plt.show()
