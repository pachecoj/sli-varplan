from multiprocessing import Pool, cpu_count
import numpy as np
import random

def do_work(map):
    i, rng = map    
    print('\tWorker: %d, Rand: %f' % (i, rng.random_sample()))
    return rng

Nworkers = cpu_count()

rngs = [np.random.RandomState(i) for i in range(Nworkers)]
pool = Pool(Nworkers)

N = 10
for n in range(N):
    print('Round: %d' % n)    
    rngs = pool.map(do_work, [(i, rngs[i]) for i in range(Nworkers)])
    
