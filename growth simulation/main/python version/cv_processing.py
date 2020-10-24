import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


def plot_cv(data, my_cmap='viridis',sigma=0, species=0):
    """
    :param data: data, 4 dims for days, parameters, sigma, species
    :param my_cmap: 'viridis', 'Greys'
    :param sigma: if dim==3, sigma = 1/10, 1/20, 1/40
    :param species: 0 for BS, 1 for No
    :return:
    """
    # attributes
    cmap = plt.get_cmap(my_cmap)
    xticks = np.linspace(0, 16, 17)
    xticklabels = ["Pf", "Nf1", "r1", "m1", "r2", "m2", "Q1c", "Q1n", "Q1p", "Q2c", "Q2n", "Q2p", "Kc1", "Kn1", "Kp1",
                  "Kn2", "Kp2"]
    yticks = np.linspace(0, 4, 5)
    yticklabels = ["5", "10", "15", "20", "25"]
    if species == 0:
        sp = 'B. Subtilis'
    else:
        sp = 'Nostoc. sp'
    #  ; sigma=5%

    fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(1, 1, 1)
    im = plt.imshow(data[:,:,sigma,species], interpolation='none', cmap=cmap)
    plt.title('Parameters\' Coefficient of Variation'+' (' + sp + ', sigma=' + str(10/np.power(2,sigma)) + '%)')
    ax.set_xlabel('Parameters')
    ax.set_ylabel('Growth Time/days')
    plt.xticks(xticks, xticklabels)  # , size=14, color='grey'
    plt.yticks(yticks, yticklabels)
    fig.colorbar(im)
    plt.show()
    # plt.savefig('cv.png')


#%% read data
mat_path = "./python version/cv_data.mat"
# mat_path = "cv_data.mat"
data = sio.loadmat(mat_path)['cv_data']


#%% plot
plot_cv(data, sigma=1, species=0)

# data[:,:,1,1]

