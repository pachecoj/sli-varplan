import sys
import numpy as np
import matplotlib.pyplot as plt
from llda_model import *
import utils
import pickle

# set random seed
np.random.seed(0)
eps = np.spacing(1)

# settings
resdir = '../data/'
name = 'SIMLLDA'
alpha = 1     # Document-level proportion Dirichlet param.
beta = 1      # Topic Dirichlet parameter
gamma = 0.1   # Label Dirichlet hyperparameter
D = 50        # Number of documents
Nd = 25       # Number of words
K = 10
W = 25
plabel = eps  # Probability of 'incorrect' annotation

# set label distribution to prefer true topic
Nl = K          # Annotation labels = Topics (don't change)
ppi = plabel * np.ones( (K, Nl) ) + np.eye(K)
ppi = ppi * (1 / np.sum(ppi, axis=1))

# init document figure
fig_doc, axs_doc = plt.subplots(5, 5)
fig_doc.suptitle('Documents')
axs_doc = axs_doc.reshape(5*5)

# draw topics and topic proportions
phi = np.zeros((K,W))
for k in range(K):
    phi[k,:] = np.random.dirichlet(beta * np.ones(W))
theta = np.random.dirichlet(alpha * np.ones(K), size=D )

# draw documents
wordcount = np.zeros((W,D))
labels = np.zeros((Nd,D), dtype='int')
z = np.zeros((D,Nd), dtype='int')
for d in range(D):

    # sample topic labels
    z_d = np.random.choice(K, (Nd,), p=theta[d,:])
    z[d,:] = z_d

    # sample words / labels
    w_d = np.zeros((Nd,), dtype='int')
    l_d = np.zeros((Nd,), dtype='int')
    for n in range(Nd):
        w_d[n] = np.random.choice(W, (1,), p=phi[z_d[n],:])
        l_d[n] = np.random.choice(Nl, (1,), p=ppi[z_d[n],:])

    # compute wordcounts
    this_wordcount, tmp = np.histogram(w_d, range(W+1))
    wordcount[:,d] = this_wordcount

    # aggregate labels
    idx_sorted = np.argsort(w_d)
    labels[:,d] = l_d[idx_sorted]

    # show document
    if d<25:
        axs_doc[d].imshow(
            np.reshape( this_wordcount, (5,5) ), cmap="gray" )
        axs_doc[d].set_xticks([])
        axs_doc[d].set_yticks([])

# show topics
utils.show_bars_topics(phi, W, K)
plt.show()

# write model
model = LLDA_Model(alpha, beta, gamma, K, Nl, Nd, ppi, phi, theta, z)
fname_model = resdir + name + '_model.dat'
f = open(fname_model, 'wb')
pickle.dump(model,f)
f.close()
print('Saved model: %s' % fname_model)

# write vocab / data
utils.write_vocab(resdir, name, W)
utils.write_data(resdir, name, wordcount, labels)
