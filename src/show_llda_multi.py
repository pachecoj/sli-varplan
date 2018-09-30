import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
import utils

datadir = '../data/'

files = ('corpus_SPARSEBARS_alg_llda_gibbs_sel_random_numinit_150_trainsteps_100_Nsamp_1000_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_150_trainsteps_100_Nsamp_1000_DISCVAR',         
         'corpus_SPARSEBARS_alg_llda_ep_sel_random_numinit_150_trainsteps_100_Nsamp_1000_DISCVAR',
         #'corpus_SPARSEBARS_alg_llda_ep_sel_empirical_numinit_150_trainsteps_100_Nsamp_1000_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_ep_sel_discvar_numinit_150_trainsteps_100_Nsamp_100_DEBUG_DISCVAR_SMALLSAMP',
         )
names = ('Inf: Gibbs, Plan: (Random, N=1k)',
         'Inf: Gibbs, Plan: Empirical MI (N=1k)',
         'Inf: Variational, Plan: Random',
         #'Inf: Variational, Plan: Empirical MI (N=1k)',         
         'Inf: Variational, Plan: MI (N=100)',
)
colors = ('--b',
          '-b',
          '--k',
          #'-+k',
          '-k',
)
facecolors = ('blue',
              'blue',
              'gray',
              #'gray',
              'gray',
)
          
Nruns = 10

# plot results
fig_tv = plt.figure()
fig_ent = plt.figure()
for fname, algname, idx in zip(files, names, range(len(files))):
    for runid in range(Nruns):    
        this_fname = '%s%s_runid_%d.mat' % (datadir, fname, runid)
        print(this_fname)
        res = sio.loadmat( this_fname )

        # concatenate errors
        if runid==0:
            errs = res['err'].flatten()[:,np.newaxis]
        else:
            errs = np.concatenate((errs,res['err'].flatten()[:,np.newaxis]), axis=1)

        # plot entropy
        plt.figure( fig_ent.number )
        Hreal = res['Hreal'].flatten()
        if (runid==0):
            plt.plot( range(len(Hreal)), Hreal, colors[idx], label=algname)
        else:
            plt.plot( range(len(Hreal)), Hreal, colors[idx])
            
    # plot TV error
    plt.figure( fig_tv.number )
    err_mean = np.mean(errs, axis=1)
    err_std = np.std(errs, axis=1)
    plt.plot( range(len(err_mean)), err_mean, colors[idx], label=algname)
    plt.fill_between( range(len(err_mean)), err_mean-err_std, err_mean+err_std, facecolor=facecolors[idx], alpha='0.3')

# decorate figures
plt.figure( fig_tv.number )
plt.legend()
plt.xlabel('Training Epoch')
plt.ylabel('Topic Error (TV)')
plt.figure( fig_ent.number )
plt.legend()
plt.xlabel('Training Epoch')
plt.ylabel('Posterior Topic Entropy')
plt.show()
