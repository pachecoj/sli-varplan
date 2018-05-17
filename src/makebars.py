import sys
import numpy as np
import matplotlib.pyplot as plt
from llda_model import *
import utils
import pickle

def write_vocab(resdir, name, V):
    fname = resdir + name + '_vocab.dat'
    f = open(fname, 'w')
    for v in range(V):
        f.write('word_' + repr(v) + '\n')
    f.close()
    print('Saved vocab: %s' % fname)

def write_data(resdir, name, wordcount, labels):
    V, D = np.shape( wordcount )
    Nd = np.shape( labels )
    Nd = Nd[0]

    # write words
    fname_words = resdir + name + '_train.dat'
    f = open(fname_words, 'w')
    for d in range(D):
        f.write( "%d" % np.sum( wordcount[:,d] ) )
        for w in range(V):
            if wordcount[w,d]>0:
                f.write( " %d:%d" % (w, wordcount[w,d]) )
        f.write('\n')
    f.close()
    print('Saved words: %s' % fname_words)

    # write labels
    fname_labels = resdir + name + '_labels.dat'
    f = open(fname_labels, 'w')
    for d in range(D):
        for n in range(Nd):
            f.write("%d " % labels[n,d])
        f.write('\n')
    f.close()
    print('Saved labels: %s' % fname_labels)

# set random seed
np.random.seed(0)
eps = np.spacing(1)

# settings
resdir = '../data/'
name = 'TINYBARS'
alpha = 1     # Document-level proportion Dirichlet param.
beta = 1      # Topic Dirichlet parameter
gamma = 0.1   # Label Dirichlet hyperparameter
D = 50        # Number of documents
Nd = 25       # Number of words
plabel = eps  # Probability of 'incorrect' annotation

# define "bars" topic dimensions
bar_width = 1    # Bar width
num_bars = 5     # Number of bars in each topic
K = 2*num_bars   # Number of topics (don't change)
V = (bar_width*num_bars)**2   # Vocabulary size (don't change)

# "vertical" topics
phi = np.zeros((K,V))
start_col = 0
for n in range(num_bars):
    this_topic = eps + np.zeros( (bar_width*num_bars, bar_width*num_bars) )
    end_col = start_col + bar_width
    this_topic[ :, start_col:end_col ] = 1
    phi[n,:] = this_topic.flatten() / np.sum( this_topic.flatten() )
    start_col = end_col

# "horizontal" topics
start_row = 0
for n in range(num_bars):
    this_topic = eps + np.zeros( (bar_width*num_bars, bar_width*num_bars) )
    end_row = start_row + bar_width
    this_topic[ start_row:end_row, : ] = 1
    phi[(num_bars+n),:] = this_topic.flatten() / np.sum( this_topic.flatten() )
    start_row = end_row

# document-level topic proportions
theta = np.random.dirichlet( alpha * np.ones((1,K)).flatten(), D )

# set label distribution to prefer true topic
Nl = K          # Annotation labels = Topics (don't change)
ppi = plabel * np.ones( (K, Nl) ) + np.eye(K)
ppi = ppi * (1 / np.sum(ppi, axis=1))

# init document figure
fig_doc, axs_doc = plt.subplots(5, 5)
fig_doc.suptitle('Documents')
axs_doc = axs_doc.reshape(5*5)

# sample documents
wordcount = np.zeros((V,D));  labels = np.zeros((Nd,D), dtype='int')
z = np.zeros((D,Nd), dtype='int')
for d in range(D):

    # sample topic labels
    z_d = np.random.choice(K, (Nd,), p=theta[d,:])
    z[d,:] = z_d

    # sample words / labels
    w_d = np.zeros((Nd,), dtype='int')
    l_d = np.zeros((Nd,), dtype='int')
    for n in range(Nd):
        w_d[n] = np.random.choice(V, (1,), p=phi[z_d[n],:])
        l_d[n] = np.random.choice(Nl, (1,), p=ppi[z_d[n],:])

    # compute wordcounts
    this_wordcount, tmp = np.histogram(w_d, range(V+1))
    wordcount[:,d] = this_wordcount

    # aggregate labels
    idx_sorted = np.argsort(w_d)
    labels[:,d] = l_d[idx_sorted]

    # show document
    if d<25:
        axs_doc[d].imshow(
            np.reshape( this_wordcount, [bar_width*num_bars, bar_width*num_bars] ), cmap="gray" )
        axs_doc[d].set_xticks([])
        axs_doc[d].set_yticks([])


# show topics
utils.show_bars_topics(phi, V, K)
plt.show()

# write model
model = LLDA_Model(alpha, beta, gamma, K, Nl, Nd, ppi, phi, theta, z)
fname_model = resdir + name + '_model.dat'
f = open(fname_model, 'wb')
pickle.dump(model,f)
f.close()
print('Saved model: %s' % fname_model)

# write vocab / data
write_vocab(resdir, name, V)
write_data(resdir, name, wordcount, labels)
