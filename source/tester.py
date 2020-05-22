import time
import multiprocessing
import random
import numpy as np
def basic_func(x):
    if x == 0:
        return 'zero'
    elif x%2 == 0:
        return 'even'
    else:
        return 'odd'

def multiprocessing_func(x, y):
    z = x*y
    time.sleep(y)
    print('{} times {} results in a/an {} number'.format(x, y, basic_func(z)))
    return z, x

if __name__ == '__main__':

    starttime = time.time()
    Ncores = multiprocessing.cpu_count() - 2
    pool = multiprocessing.Pool(Ncores)
    list = []
    for i in range(Ncores):
        list.append((i, np.random.randint(1,Ncores)))
    print(list)



    print('Returns:')
    print(pool.starmap(multiprocessing_func, list))
    pool.close()
    print('That took {} seconds'.format(time.time() - starttime))
