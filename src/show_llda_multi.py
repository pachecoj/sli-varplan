import matplotlib.pyplot as plt
import scipy.io as sio
import utils
import ipdb

datadir = '../data/'

# files = ('corpus_TINYBARS_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_TINYBARS_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_Nsamp_1000_DBG',
#          'corpus_TINYBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DBG')
# names = ('Empirical (N=1K)','Variational','Random')
# colors = ('-r','-k','-b')
# Nruns = 10

files = ('corpus_TINYBARS_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_Nsamp_1000_AVG',
         'corpus_TINYBARS_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_Nsamp_1000',
         'corpus_TINYBARS_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_DBG',
         'corpus_TINYBARS_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_Nsamp_1000_DBG',
         'corpus_TINYBARS_alg_llda_ep_sel_hybrid_numinit_50_trainsteps_100_Nsamp_1000_DBG',
         'corpus_TINYBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_DBG')
names = ('Gibbs (Mean, N=1K)','Gibbs (MAP, N=1K)','Empirical (N=1K)','Variational','Hybrid (N=1K)','Random')
colors = ('-g','--g','-r','-k','--k','-b')
Nruns = 10

# plot results
fig_tv = plt.figure()
for runid in range(Nruns):
    for fname, algname, idx in zip(files, names, range(len(files))):
        this_fname = '%s%s_runid_%d.mat' % (datadir, fname, runid)
        print(this_fname)
        res = sio.loadmat( this_fname )
        err = res['err'].flatten()
        plt.figure( fig_tv.number )
        if (runid == 0):
            plt.plot( range(len(err)), err, colors[idx], label=algname)
        else:
            plt.plot( range(len(err)), err, colors[idx])


plt.figure( fig_tv.number )
plt.legend()
plt.xlabel('Training Epoch')
plt.ylabel('Topic Error (TV)')
plt.show()
