# try CA version 2, Modelling the morphology of migrating bacterial colonies
# we use rectangular grid or hexagonal grid
# dtype='uint16', 0~65535

from utils3 import *
from initial_values import *

global ELE
dtype = 'int16'

# %% simulation
# fundamental parameters
'''elements of a state (6 by n by n)
0 for b of BS
1 for s of BS
2 for b of No
3 for s of No
4 for nutrition
5 for EPS content
'''
n = 45  # size
# n = 70
# epoch = 130  # number of period. more causes NaN??
epoch = 100

# probs_migration, weight, params_BS, params_No, init_classic = v0_10_22()
# probs_migration, weight, params_BS, params_No, init_classic = v1_high_no()
# probs_migration, weight, params_BS, params_No, init_classic = v2_low_nut()
# probs_migration, weight, params_BS, params_No, init_classic = v4_high_rps()

probs_migration, weight, params_BS, params_No, init_classic = v3_to_nut()
state_init = init_directed_nutrient(n, num=init_classic, dtype=dtype)

# probs_migration, weight, params_BS, params_No, init_classic = v2_low_nut()
# state_init = init_random_nutrient(n, num=init_classic, dtype=dtype)

# n1 = int(n/3)
# n2 = int(n/3*2)
# state_init = init_multi_pos(n, [(n1,n1),(n1,n2),(n2,n2),(n2,n1)], np.repeat(np.array([[100,200]],dtype=dtype), 4, axis=0), 200, 1000)


# alpha1=1.4071;
# beta1=0.0093;
# alpha2=10.1682;
# beta2=0.0057;
# the effect: just grow?
'''
notes:
- EPS has a far more larger number (a unit for less...); weight contain info about unit 
    unificaiton, so we set a low value
'''

# the init_simple2d added a dimension to a state, which fits the original evolution function
# put a state 2 cell in the center, but we can write our own
# each time I run the simulation, state_init changed ????! have to put it here
# state_init = init_classic_center(n, num=init_classic, dtype=dtype)
states_epoch, states_phase = stimulation_v3(n=n, state_init=state_init, grid='rect_Moore',
                                            # rect_Moore, rect_Neumann, hexagonal
                                            epoch=epoch, see_phase=False,
                                            probs_migration=probs_migration, weight=weight,
                                            params_BS=params_BS, params_No=params_No,
                                            dtype=dtype)
# it is still hard to perform a large-scale simulation


#%%
def my_plot2d_animate(ca, idx=0, title='evolve', interval=50, my_cmap=None):
    fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(1, 1, 1)
    font_title = {'family': 'Arial',  # 'Arial', 'times new roman', 'sans-serif', 'Cambria Math', 'Comic Sans MS'
                  'weight': 1.5,
                  'size': 16, }
    font_tick = {'family': 'Arial',  # 'Arial', 'times new roman', 'sans-serif', 'Cambria Math', 'Comic Sans MS'
                  'weight': 'normal',
                  'size': 12, }
    # cmap setting
    if my_cmap is None:
        my_cmap = 'viridis'  # default
    if my_cmap == 'bg':  # brown to green
        # colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # rgb
        colors = [(1, 204 / 255, 0), (0, 1, 0)]  # standard brown and green
        # colors = [(239 / 255, 224 / 255, 203 / 255), (214 / 255, 242 / 255, 217 / 255)]  # get from our wiki
        # colors = [(239 / 255, 224 / 255, 203 / 255), (214 / 255, 242 / 255, 217 / 255)]  # try to adjust
        my_cmap = LinearSegmentedColormap.from_list(name='bg', colors=colors)  # , N=n_bin
    cmap = plt.get_cmap(my_cmap)  # color
    fig.set_facecolor((218/255, 238/255, 214/255))

    im = plt.imshow(ca[0][idx], animated=True, cmap=cmap, vmin=0, vmax=np.max(ca[:, idx]))  # set the value range here!
    fig.colorbar(im)  # , ax=ax; using default color bar setting, but will not change with time to con
    # unsuccessful, our standard ticks, set special values in the bar
    # plt.colorbar(im, cmap=my_cmap, ticks=get_ticks(np.max(ca[:, idx])))
    a = plt.gca()
    xlabel = a.get_xticklabels(im)
    ylabel = a.get_yticklabels(im)
    a.set_xticklabels(xlabel, fontdict=font_tick)
    a.set_yticklabels(ylabel, fontdict=font_tick)
    i = {'index': 0}  # because it is used in a newly-defined function

    def updatefig(*args):
        i['index'] += 1
        if i['index'] == len(ca):
            i['index'] = 0
        title = 'evolve at epoch ' + str(i['index']) + ' for index ' + str(idx)  # manually set 'epoch' or 'stage'
        ax.set_title(title, fontdict=font_title)
        im.set_array(ca[i['index']][idx])
        # plt.colorbar(cmap=my_cmap, ticks=get_ticks(np.max(ca[:, idx])))
        # fig.colorbar(im)
        # plt.rcParams['savefig.facecolor'] = (215 / 255, 239 / 255, 214 / 255)
        plt.rcParams['savefig.facecolor'] = (251 / 255, 238 / 255, 201 / 255)
        return im,

    ani = animation.FuncAnimation(fig, updatefig, interval=interval, blit=True)
    plt.show()
    # save as a gif file, changing writer resolved the error
    # have to enlarge maximum memory limit
    ani.save(filename="evolve.gif", writer='pillow')  # ,dpi=60


# %% result
# a dynamic graph, repeat playing
idx = 0
# each epoch costs "interval" millisecond
# my_plot2d_animate(states_epoch, idx=idx, interval=200, my_cmap='Greys')
my_plot2d_animate(states_epoch, idx=idx, interval=100)  # None=='viridis', **high contrast**
# my_plot2d_animate(states_epoch, idx=idx, interval=200, my_cmap='bg')  # brown green, project
# other choices: mpl.cm.bwr,


# %% for debugging purposes
time = 66
# time = -1
idx = 4
my_plot2d(states_epoch, timestep=time, idx=idx)
# states_epoch[time]
# states_epoch[time][idx]
# my_plot2d(state_update, idx=0)
# my_plot2d_animate(states_phase)  # dynamic graph by phase
# my_plot2d(states_phase, timestep=2)
# my_plot2d(state_init, idx=0)
# my_plot2d(num_rec2)
# my_plot2d(num_spo1)
# my_plot2d(grow)

