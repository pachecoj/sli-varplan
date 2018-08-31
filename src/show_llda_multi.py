import matplotlib.pyplot as plt
import scipy.io as sio
import utils

datadir = '../data/'

files = ('corpus_SPARSEBARS_alg_llda_gibbs_sel_random_numinit_150_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_150_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_ep_sel_random_numinit_150_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_ep_sel_empirical_numinit_150_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',
         # 'corpus_SPARSEBARS_alg_llda_ep_sel_discvar_numinit_150_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR2',
         'corpus_SPARSEBARS_alg_llda_ep_sel_discvar_numinit_150_trainsteps_100_Nsamp_100_DEBUG_DISCVAR_SMALLSAMP',               
         )
names = ('Inf: Gibbs, Plan: Random)',
         'Inf: Gibbs, Plan: Empirical MI (N=1k)',
         'Inf: Variational, Plan: Random',
         'Inf: Variational, Plan: Empirical MI (N=1k)',         
         'Inf: Variational, Plan: MI (N=1k, Discrim.)',
)
colors = ('--b',
          '-+b',
          '--k',
          '-+k',
          '-k',
)
Nruns = 1

# plot results
fig_tv = plt.figure()
fig_ent = plt.figure()
for runid in range(Nruns):
    for fname, algname, idx in zip(files, names, range(len(files))):
        this_fname = '%s%s_runid_%d.mat' % (datadir, fname, runid)
        print(this_fname)
        res = sio.loadmat( this_fname )

        # plot TV error
        plt.figure( fig_tv.number )
        err = res['err'].flatten()
        if (runid == 0):
            plt.plot( range(len(err)), err, colors[idx], label=algname)
        else:
            plt.plot( range(len(err)), err, colors[idx])

        # plot entropy
        plt.figure( fig_ent.number )
        Hreal = res['Hreal'].flatten()
        if (runid==0):
            plt.plot( range(len(Hreal)), Hreal, colors[idx], label=algname)
        else:
            plt.plot( range(len(Hreal)), Hreal, colors[idx])


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
