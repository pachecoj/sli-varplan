import matplotlib.pyplot as plt
import scipy.io as sio
import utils
import ipdb

datadir = '../data/'

files = ('corpus_SPARSEBARS_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_50_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',
         'corpus_SPARSEBARS_alg_llda_ep_sel_discvar_numinit_50_trainsteps_100_Nsamp_1000_DEBUG_DISCVAR',         
         )
names = ('Inf: Gibbs, Plan: Random)',
         'Inf: Gibbs, Plan: MI (N=1k)',
         'Inf: Variational, Plan: Random',
         'Inf: Variational, Plan: MI (N=1k, Discrim.)',
)
colors = ('--b',
          '-b',
          '--k',
          '-k',
)
Nruns = 2

# files = ('corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_50_trainsteps_100_Nsamp_1000_FIXENTROPY',
#          'corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_50_trainsteps_100_Nsamp_2000_HIGHER_SAMPBURN',
#          'corpus_SPARSEBARS_alg_llda_gibbs_sel_gibbs_numinit_50_trainsteps_100_Nsamp_4000_HIGHER_SAMPBURN2',
#          'corpus_SPARSEBARS_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_SPARSEBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_SPARSEBARS_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_DBG'
#          )
# names = ('Gibbs (N=1k, Burn=100)',
#     'Gibbs (N=2k, Burn=200)',
#     'Gibbs (N=4k, Burn=200)',
#     'Variational',
#     'Random',
#     'Empirical (N=1K)'
#     )
# colors = ('-g',
#     '--g',
#     '-+g',
#     '-k',
#     '-b',
#     '-r'
#     )
# Nruns = 1

# files = ('corpus_SIMLLDA_alg_llda_gibbs_sel_gibbs_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_SIMLLDA_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DBG_ALIGNTRUE',
#          'corpus_SIMLLDA_alg_llda_gibbs_sel_gibbs_numinit_50_trainsteps_100_Nsamp_1000_FIXENTROPY',
#          'corpus_SIMLLDA_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_SIMLLDA_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_SIMLLDA_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_DBG'
#          )
# names = ('Gibbs (Align MAP, MI)',
#     'Gibbs (Align GT, Random)',
#     'Gibbs (Align GT, MI)',
#     'Variational',
#     'Random',
#     'Empirical (N=1K)'
#     )
# colors = ('-g',
#     '--g',
#     '-+g',
#     '-k',
#     '-b',
#     '-r'
#     )
# Nruns = 3

#          'corpus_TINYBARS_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_Nsamp_1000_REALENT',
#          'corpus_TINYBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_REALENT',
#          'corpus_TINYBARS_alg_llda_ep_sel_hybrid_numinit_50_trainsteps_100_Nsamp_1000_REALENT',
#          'corpus_TINYBARS_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_REALENT')
# names = ('Gibbs (N=1K)','Variational','Random','Hybrid (N=1K)','Empirical (N=1K)')
# colors = ('-g','-k','-b','--k','-r')
# Nruns = 10

# files = ('corpus_TINYBARS_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_TINYBARS_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_TINYBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DBG')
# names = ('Empirical (N=1K)','Variational','Random')
# colors = ('-r','-k','-b')
# Nruns = 10

# files = ('corpus_TINYBARS_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_Nsamp_1000_AVG',
#          'corpus_TINYBARS_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_Nsamp_1000',
#          'corpus_TINYBARS_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_TINYBARS_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_TINYBARS_alg_llda_ep_sel_hybrid_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_TINYBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DBG')
# names = ('Gibbs (Mean, N=1K)','Gibbs (MAP, N=1K)','Empirical (N=1K)','Variational','Hybrid (N=1K)','Random')
# colors = ('-g','--g','-r','-k','--k','-b')
# Nruns = 10

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
