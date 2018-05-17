import matplotlib.pyplot as plt
import scipy.io as sio
import utils

datadir = '../data/'

# files = ('corpus_BARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_runid_0.mat',
#     'corpus_BARS_alg_llda_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_1000_runid_0.mat',
#     'corpus_BARS_alg_llda_ep_sel_variational_numinit_50_trainsteps_100_runid_0.mat')
# names = ('Random', 'Empirical (N=1k)', 'Variational')
# colors = ('-b','-r','-k')

# files = ('corpus_TINYBARS_alg_llda_batch_ep_sel_random_numinit_50_trainsteps_100_runid_0.mat',
#          'corpus_TINYBARS_alg_llda_batch_ep_sel_empirical_numinit_50_trainsteps_100_Nsamp_5000_runid_0.mat',
#          'corpus_TINYBARS_alg_llda_batch_ep_sel_empiricalindependent_numinit_50_trainsteps_100_Nsamp_5000_runid_0.mat',
#     'corpus_TINYBARS_alg_llda_batch_ep_sel_hybrid_numinit_50_trainsteps_100_Nsamp_5000_runid_0.mat',
#     'corpus_TINYBARS_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_runid_0.mat')
# names = ('Random', 'Empirical (N=5k, Dep.)', 'Empirical (N=5k, Indep.)', 'Hybrid (N=5k)','MCMC (N=1k)')
# colors = ('-b','-r','--r','-k','-g')

files = ('corpus_TINYBARS_alg_llda_ep_sel_random_numinit_50_trainsteps_100_Nsamp_1000_runid_0.mat',         
    'corpus_TINYBARS_alg_llda_gibbs_sel_random_numinit_50_trainsteps_100_Nsamp_1000_runid_0.mat')
names = ('Random','MCMC (N=1k)')
colors = ('-b','-g')

# plot results
fig_tv = plt.figure()
for fname, algname, idx in zip(files, names, range(len(files))):
    print(datadir + fname)
    res = sio.loadmat( datadir + fname )
    err = res['err'].flatten()
    plt.figure( fig_tv.number )
    plt.plot( range(len(err)), err, colors[idx], label=algname)

    # plot topics
    this_lambda = res['lambda']
    fig_top, ax_top = utils.show_bars_topics(this_lambda[:,:,-1], int(res['W'][0]), int(res['K'][0]))
    fig_top.suptitle('%s Topics' % algname)


plt.figure( fig_tv.number )
plt.legend()
plt.xlabel('Training Epoch')
plt.ylabel('Topic Error (TV)')
plt.show()
