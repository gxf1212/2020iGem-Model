import numpy as np


# define useful functions
def MM(r, k, n):
    return r / (r + k * n)


def toxin(ga, n1, n2):
    return 1 - ga * n1 * n2


def save_files(N1, N2, Rc, Rn, Rp, time, G1, G2):
    np.savetxt('./main/data/N1.csv', N1, delimiter=',')
    np.savetxt('./main/data/N2.csv', N2, delimiter=',')
    np.savetxt('./main/data/Rc.csv', Rc, delimiter=',')
    np.savetxt('./main/data/Rn.csv', Rn, delimiter=',')
    np.savetxt('./main/data/Rp.csv', Rp, delimiter=',')
    np.savetxt('./main/data/time.csv', time, delimiter=',')
    np.savetxt('./main/data/G1.csv', G1, delimiter=',')
    np.savetxt('./main/data/G2.csv', G2, delimiter=',')



