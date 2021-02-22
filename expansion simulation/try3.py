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
# n = 25  # size
n = 40
epoch = 130  # number of period. more causes NaN??
# epoch = 20

probs_migration, weight, params_BS, params_No, init_classic = v0_10_22()
# probs_migration, weight, params_BS, params_No, init_classic = v1_high_no()
# probs_migration, weight, params_BS, params_No, init_classic = v2_low_nut()

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
state_init = init_classic_center(n, num=init_classic, dtype=dtype)
# n1 = int(n/3)
# n2 = int(n/3*2)
# state_init = init_multi_pos(n, [(n1,n1),(n1,n2),(n2,n2),(n2,n1)], np.repeat(np.array([[100,200]],dtype=dtype), 4, axis=0), 200, 1000)

# probs_migration, weight, params_BS, params_No, init_classic = v2_low_nut()
# state_init = init_random_nutrient(n, num=init_classic, dtype=dtype)

states_epoch, states_phase = stimulation_v3(n=n, state_init=state_init, grid='rect_Moore',
                                            # rect_Moore, rect_Neumann, hexagonal
                                            epoch=epoch, see_phase=False,
                                            probs_migration=probs_migration, weight=weight,
                                            params_BS=params_BS, params_No=params_No,
                                            dtype=dtype)
# it is still hard to perform a large-scale simulation

# %% load data
states_epoch = np.load('multipoint.npy')

# %% result
# a dynamic graph, repeat playing
idx = 0
# each epoch costs "interval" millisecond
# my_plot2d_animate(states_epoch, idx=idx, interval=200, my_cmap='Greys')
my_plot2d_animate(states_epoch, idx=idx, interval=50)  # None=='viridis', **high contrast**
# my_plot2d_animate(states_epoch, idx=idx, interval=200, my_cmap='bg')  # brown green, project
# other choices: mpl.cm.bwr,

#%% set a range to write a gif
idx = 0
# states = states_epoch[:,:,5:45,5:45]
# my_plot2d_animate(states, idx=idx, interval=50)
# %% for debugging purposes
time = 46
# time = -1
idx = 0
my_plot2d(states_epoch, timestep=time, idx=idx, my_cmap='viridis')
# states_epoch[time]
# states_epoch[time][idx]
# my_plot2d(state_update, idx=0)
# my_plot2d_animate(states_phase)  # dynamic graph by phase
# my_plot2d(states_phase, timestep=2)
# my_plot2d(state_init, idx=0)
# my_plot2d(num_rec2)
# my_plot2d(num_spo1)
# my_plot2d(grow)

