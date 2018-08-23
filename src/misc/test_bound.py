import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
from scipy.optimize import minimize, fmin_l_bfgs_b
import scipy.stats as stats

np.random.seed(100000001)

def unwrap(x,K):
    eta = np.array(x[0:K])
    lam = np.array(x[K:])
    return eta, lam

def f(x, *args):
    K = args[0]
    y_samp = args[1]
    x_samp = args[2]
    N = y_samp.shape[0]
    eta, lam = unwrap(x,K)
    v = lam**(-1)
    m = v * eta

    # import ipdb
    # ipdb.set_trace()

    obj = 0.
    for n in range(N):
        k = y_samp[n]
        this_x = x_samp[n]
        obj -= 1/N * stats.norm.logpdf(this_x, loc=m[k], scale=np.sqrt(v[k]))

    return obj

def J(x, *args):
    K = args[0]
    y_samp = args[1]
    x_samp = args[2]
    N = y_samp.shape[0]
    eta, lam = unwrap(x,K)
    v = lam**(-1)
    m = v * eta

    # compute moments
    Ex_p = np.zeros(K)
    Exx_p = np.zeros(K)
    Ex_q = np.zeros(K)
    Exx_q = np.zeros(K)
    for k in range(K):
        idx = np.where(y_samp == k )
        this_x = x_samp[idx]

        # compute model moments
        Ex_p[k] = np.mean( this_x )
        Exx_p[k] = 1/2 * np.var( this_x ) + 1/2 * Ex_p[k]**2

        # expected moments under q
        N_k = this_x.shape[0]
        Ex_q[k] = N_k / N * m[k]
        Exx_q[k] = N_k / N * 1/2 * ( v[k] + m[k]**2 )

    Jac_p = np.concatenate((Ex_p, Exx_p), axis=0)
    Jac_q = np.concatenate((Ex_q, Exx_q), axis=0)

    return Jac_p - Jac_q

# constants
eps = 1e-6

# Define GMM
# y ~ Cat(.)
# x | y ~ N(m_y, s)
K=5
p_y = np.random.random_sample(size=(K,))
p_y = p_y / np.sum( p_y )
m = (-5,-1,0,1,5)
s = (0.5, 1.0, 1.5, 2.0, 2.5)

# draw joint samples
N = 1000
y_samp = np.zeros(N, dtype='int')
x_samp = np.zeros(N)
for n in range(N):
    y_samp[n] = np.random.choice(K, p=p_y)
    x_samp[n] = np.random.normal(loc=m[ y_samp[n] ], scale=s[ y_samp[n] ])

#res = minimize(f, np.ones(2*K,), method='BFGS', jac=J, args=(K,y_samp,x_samp), options={'disp': True})
bounds = [(None,None),(None,None),(None,None),(None,None),(None,None),(eps,None),(eps,None),(eps,None),(eps,None),(eps,None)]
x, f_val, d = fmin_l_bfgs_b(f, np.ones(2*K), fprime=J, args=(K,y_samp,x_samp), approx_grad=False, bounds=bounds, iprint=1)

# print results
eta, lam = unwrap(x,K)
v = lam**(-1)
m = v * eta
print( m )
print( np.sqrt(v) )
