import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


def plot_cv(data, my_cmap='viridis', sigma=0, species=0, high=5, ticks=None, ignore=True, idx=None):
    """
    :param data: data, 4 dims for days, parameters, sigma, species
    :param my_cmap: 'viridis', 'Greys'
    :param high: percentage of sigma
    :param sigma: if dim==3, sigma = 1/10, 1/20, 1/40 (older version)
    if dim==1, sigma=1/20
    :param species: 0 for BS, 1 for No
    :param idx: the range of variables
    :return:
    """
    if idx==None:
        idx = np.linspace(0, data.shape[2], data.shape[2]+1)  # all

    # attributes
    cmap = plt.get_cmap(my_cmap)
    if ignore == True:
        data = np.delete(data, 3, axis=1)
        xticks = np.linspace(0, 15, 16)
        xticklabels = ["Pf", "Nf1", "r1", "r2", "m2", "Q1c", "Q1n", "Q1p", "Q2c", "Q2n", "Q2p", "Kc1", "Kn1",
                       "Kp1", "Kn2", "Kp2"]
    else:
        xticks = np.linspace(0, 16, 17)
        xticklabels = ["Pf", "Nf1", "r1", "m1", "r2", "m2", "Q1c", "Q1n", "Q1p", "Q2c", "Q2n", "Q2p", "Kc1", "Kn1", "Kp1",
                   "Kn2", "Kp2"]
    # xticklabels = ["$Pf$", "$Nf_1$", "$r_1$", "$m_1$", "$r_2$", "$m_2$", "$Q_{1c}$", "$Q_{1n}", "$Q_{1p}", "$Q_{2c}",
    #                "$Q_{2n}$", "$Q_{2p}$", "$K_{c1}$", "$K_{n1}$", "$K_{p1}$", "$K_{n2}$", "$K_{p2}$"]
    yticks = np.linspace(0, 4, 5)
    yticklabels = ["5", "10", "15", "20", "25"]
    if species == 0:
        sp = 'B. Subtilis'
    else:
        sp = 'Nostoc. sp'
    # title = 'Parameters\' Coefficient of Variation' + ' (' + sp + ', sigma=' + str(high / np.power(2, sigma)) + '%)'
    title = 'Parameters\' Coefficient of Variation' + ' (' + sp + ', sigma=' + str(high / np.power(2, sigma)) + '%)'
    print(title)
    # font setting
    # plt.rcParams.update({
    #     "text.usetex": True,
    #     "font.family": "Arial"})
    font_la = {'family': 'Arial',
               'weight': 'bold',
               'size': 16,
               }
    font_title = {'family': 'Arial',
                  'weight': 'bold',
                  'size': 20,
                  }
    font_tick = {'family': 'Comic Sans MS',  # 'Arial', 'times new roman', 'sans-serif', 'Cambria Math'
                 'weight': 1.5, 'size': 12}

    fig = plt.figure(figsize=(16, 5))
    im = plt.imshow(data[:, :, sigma, species], interpolation='none', cmap=cmap)
    # plt.title(title, fontdict=font_title)
    plt.xlabel('Parameters', font_la)
    plt.ylabel('Growth Time/days', font_la)
    plt.xticks(xticks, xticklabels)  # color='grey'
    plt.yticks(yticks, yticklabels)
    a = plt.gca()
    a.set_xticklabels(xticklabels, font_tick)
    a.set_yticklabels(yticklabels, font_tick)
    if ticks is not None:
        fig.colorbar(im, ticks=ticks)
    else:
        fig.colorbar(im)
    plt.show()
    return title


def plot_cv2(data, my_cmap='viridis', sigma=0, species=0, high=5, ticks=None):
    """
    :param data: data, 4 dims for days, parameters, sigma, species
    :param my_cmap: 'viridis', 'Greys'
    :param high: percentage of sigma
    :param sigma: if dim==3, sigma = 1/10, 1/20, 1/40 (older version)
    if dim==1, sigma=1/20
    :param species: 0 for BS, 1 for No
    :return:
    """
    # attributes
    cmap = plt.get_cmap(my_cmap)
    data = np.delete(data, 3, axis=1)
    xticks = np.linspace(0, 10, 11)
    xticklabels = ["Pf", "Nf1", "r2", "m1", "m2", "Q1c", "Q1n", "Q1p", "Q2c", "Q2n", "Q2p"]
    yticks = np.linspace(0, 4, 5)
    yticklabels = ["5", "10", "15", "20", "25"]
    if species == 0:
        sp = 'B. Subtilis'
    else:
        sp = 'Nostoc. sp'
    # title = 'Parameters\' Coefficient of Variation' + ' (' + sp + ', sigma=' + str(high / np.power(2, sigma)) + '%)'
    title = 'Parameters\' Coefficient of Variation' + ' (' + sp + ', sigma=' + str(high / np.power(2, sigma)) + '%)'
    print(title)
    # font setting
    # plt.rcParams.update({
    #     "text.usetex": True,
    #     "font.family": "Arial"})
    font_la = {'family': 'Arial',
               'weight': 'bold',
               'size': 16,
               }
    font_title = {'family': 'Arial',
                  'weight': 'bold',
                  'size': 20,
                  }
    font_tick = {'family': 'Comic Sans MS',  # 'Arial', 'times new roman', 'sans-serif', 'Cambria Math'
                 'weight': 1.5, 'size': 12}

    fig = plt.figure(figsize=(16, 5))
    im = plt.imshow(data[:, :, sigma, species], interpolation='none', cmap=cmap)
    # plt.title(title, fontdict=font_title)
    plt.xlabel('Parameters', font_la)
    plt.ylabel('Growth Time/days', font_la)
    plt.xticks(xticks, xticklabels)  # color='grey'
    plt.yticks(yticks, yticklabels)
    a = plt.gca()
    a.set_xticklabels(xticklabels, font_tick)
    a.set_yticklabels(yticklabels, font_tick)
    if ticks is not None:
        fig.colorbar(im, ticks=ticks)
    else:
        fig.colorbar(im)
    plt.show()
    return title


# %% read data
mat_path = "cv_data.mat"
data = sio.loadmat(mat_path)['cv_data']
mat_path = "./cv_data_v0.mat"
data0 = sio.loadmat(mat_path)['cv_data']

# %% plot
title = plot_cv(data, sigma=0, species=0, high=5, ignore=False)
# data[:,:,0,1]
# plot_cv(data0, sigma=1, species=1, high=10)

# %% final
# title = plot_cv(data, sigma=0, species=0, high=5)  # , ticks=[0, 0.2, 0.4, 0.6, 0.8, 1, 1.2]
title = plot_cv2(data[:, 0:12, :, :], sigma=0, species=0, high=5)

#%% save
plt.savefig(title+'.png', format='png', bbox_inches='tight', transparent=True, dpi=200)


