import numpy as np
from scipy import stats
from scipy.stats import multivariate_normal as mvn
import matplotlib.pyplot as plt

DEBUG=0

# def mvnlogpdf(x, m, V):
#     k = float( x.shape[1] )
#     Lam = np.linalg.inv( V )
#     D = (x - m).dot( Lam )
#     logp = -0.5 * np.sum( D * (x - m), axis=1 )
#     logp += - k / 2. * np.log( 2 * np.pi ) - 0.5 * np.log( np.linalg.det( V ) )
#     return logp

# def mvnpdf(x, m, V):
#     logp = mvnlogpdf(x, m, V)
#     return np.exp( logp )

# set random seed
np.random.seed(0)

# numerical marginal PDF
N = 1000
x = np.linspace(-100, 100, N)
dx = x[1] - x[0]
dy = dx
Nd = 50
delta = np.linspace(0, 50, Nd)
sig = 5.0

# numerical joint PDF
YY, XX = np.meshgrid(x, x)
XY = np.concatenate( ( XX.reshape((N**2,1)), YY.reshape((N**2,1)) ), 1 )

# sweep through separation values
K = np.zeros(Nd)
MI = np.zeros(Nd)
for m, i in zip( delta, range(Nd) ):
    p = 0.5 * stats.norm.pdf( x, loc=-m, scale=sig ) + 0.5 * stats.norm.pdf( x, loc=m, scale=sig )
    logp = np.log( p )

    # only use finite elements
    idx = np.isfinite( logp )
    this_x = x[ idx ]
    p = p[ idx ]
    logp = logp[ idx ]

    # moments
    mu = this_x.dot( p ) * dx
    xx = this_x * this_x
    v = xx.dot( p ) * dx - mu**2

    # projection PDF and log-PDF
    q = stats.norm.pdf( this_x, loc=mu, scale=np.sqrt( v ) )
    logq = np.log( q )

    # compute KL
    Hp = p.dot( - logp ) * dx
    Hpq = p.dot( - logq ) * dx
    K[i] = Hpq - Hp

    # compute joint PDF
    m0_xy = np.array([-m, 0])
    m1_xy = np.array([m, 0])
    Vxy = sig**2 * np.eye(2)
    pxy = 0.5 * mvn.pdf(XY, m0_xy, Vxy) + 0.5 * mvn.pdf(XY, m1_xy, Vxy)
    logpxy = np.log( pxy )

    # compute marginal entropy H(X)
    this_pxy = np.reshape( pxy, (N,N) )
    px = np.sum( this_pxy, axis=1 ) * dx
    logpx = np.log( px )
    idx = np.isfinite( logpx )
    px = px[ idx ]
    logpx = logpx[ idx ]
    Hx = px.dot( - logpx ) * dx

    # DEBUG
    if DEBUG:
        fig, axs = plt.subplots(1, 3)
        fig.set_figwidth( 17 )
        axs[0].plot( XX[idx,0], px, label='Px' )

    # compute marginal entropy H(Y)
    py = np.sum( this_pxy.T, axis=1 ) * dy
    logpy = np.log( py )
    idx = np.isfinite( logpy )
    py = py[ idx ]
    logpy = logpy[ idx ]
    Hy = py.dot( -logpy ) * dy

    # DEBUG
    if DEBUG:
        axs[0].plot( YY[0,idx], py, label='Py' )
        axs[0].legend()
        axs[0].set_xlabel('X / Y')

    # DEBUG
    if DEBUG:
        axs[1].contour(XX, YY, np.reshape(pxy, (N,N)))
        axs[1].set_xlabel('X')
        axs[1].set_ylabel('Y')

    # compute joint entropy H(X,Y)
    idx = np.isfinite( logpxy )
    pxy = pxy[ idx ]
    logpxy = logpxy[ idx ]
    Hxy = pxy.dot( -logpxy ) * dx * dy

    # compute / store MI
    MI[i] = Hx + Hy - Hxy

    # DEBUG
    if DEBUG:
        axs[2].bar( (1,2,3), (Hx, Hy, Hxy), tick_label=('Hx','Hy','Hxy') )
        plt.show()


# plot stuff
plt.plot(delta, K, label='KL(p(x)||q(x))')
plt.plot(delta, MI, label='I(X;Y)')
plt.xlabel('Mean Separation')
plt.legend()
plt.show()
