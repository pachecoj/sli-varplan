import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

plt.rc('text', usetex=True)
plt.rc('font', size=16)

x = np.linspace(-5, 5, num=1000)
dx = x[1] - x[0]
y = np.linspace(-5, 5, num=1000)

# marginal
m0 = -2.
m1 = 2. 
logp_x = np.log( 0.5 * stats.norm.pdf(x, loc=m0) + 0.5 * stats.norm.pdf(x, loc=m1) )
logp_x -= np.max(logp_x)
p_x = np.exp(logp_x) / np.sum(np.exp(logp_x)) / dx
v_q = 0.5 * (1. + m0**2) + 0.5 * (1. + m1**2) - 0.5 * m0**2 - 0.5 * m0**2
logq_x = stats.norm.logpdf(x, loc=0., scale=np.sqrt(v_q))
q_x = np.exp(logq_x)
fig = plt.figure()
plt.plot(x, p_x, '-b', label='$\displaystyle p(x\mid \mathcal{Y})$')
plt.plot(x, q_x, '-k', label='$\displaystyle q(x)$')
plt.xlabel('X')
plt.title('Posterior Inference')
plt.legend()
fig.axes[0].set_xticks([])
fig.axes[0].set_yticks([])

# compute joint
a = 2.
sinx = np.sin(a*x)
logp_y = -0.5 * (y[:,np.newaxis] - sinx[np.newaxis,:])**2
logp = logp_y + logq_x[np.newaxis,:]
logp = logp - np.max(logp)
p = np.exp(logp) / np.sum(np.exp(logp)) / dx / dx
fig = plt.figure()
plt.contour(p,10)
plt.xlabel('X')
plt.ylabel('Y')
fig.axes[0].set_xticks([])
fig.axes[0].set_yticks([])
plt.title('Local Approximation $\displaystyle\hat{p}(x,y)$')

# compute linar conditional approximation
c = 0.15
v = 0.25
p_cond = p / np.sum(p, axis=1)[:,np.newaxis] / dx
logq = -0.5 * (x[np.newaxis,:] - c*y[:,np.newaxis])**2/v
logq = logq - np.max(logq)
q = np.exp(logq) / np.sum(np.exp(logq)) / dx
fig = plt.figure()
plt.contour(p_cond,10)
plt.contour(q,10,cmap='gray')
fig.axes[0].set_xticks([])
fig.axes[0].set_yticks([])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Linear Gaussian Conditional Approximation')

# compute nonlinear conditional approximation
b = 1.5
b0 = -0.8
v2 = .5
fy = 1 / (1 + np.exp(-y))
logq2 = -0.5 * (x[np.newaxis,:] - b*fy[:,np.newaxis] - b0)**2/v2
logq2 -= np.max(logq2)
q2 = np.exp(logq2) / np.sum(np.exp(logq2)) / dx
fig = plt.figure()
plt.contour(p_cond,10)
plt.contour(q2,10,cmap='gray')
fig.axes[0].set_xticks([])
fig.axes[0].set_yticks([])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Nonlinear Gaussian Conditional Approximation')

plt.show()
