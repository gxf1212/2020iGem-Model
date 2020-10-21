# try CA version 2, Modelling the morphology of migrating bacterial colonies
# we use rectangular grid or hexagonal grid
# dtype='uint16', 0~65535

from utils3 import *

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
n = 25  # size
epoch = 100  # number of period. more causes NaN??
# epoch = 20
probs_migration = [0.02, 0.02, 0.002, 0.02]  # diffusion coefficient
weight = [[1, 0.05, 1], [1, 0.05, 1]]  # cell, EPS and nutrient on migration
# stable params
# growth, to spores, nutrient consumption, rps production (with N), carrying capacity
# with dN/dt?
# params_BS = (0.3, 0.03, 0.1, 0.5, 1000)
# params_No = (0.3, 0.01, -0.1, 0.5, 1000)

# real params
# these params should be proportional to the original
# growth/decay rate: *20
params_BS = (0.40, 0.057, 0.1, 0.5, 1000)
params_No = (0.01, 0.001, -0.1, 0.5, 1000)

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
state_init = init_classic_center(n, num=(200, 200, 100, 1000), dtype=dtype)
states_epoch, states_phase = stimulation_v3(n=n, state_init=state_init, grid='rect_Moore',
                                            # rect_Moore, rect_Neumann, hexagonal
                                            epoch=epoch, see_phase=False,
                                            probs_migration=probs_migration, weight=weight,
                                            params_BS=params_BS, params_No=params_No,
                                            dtype=dtype)
# it is still hard to perform a large-scale simulation

# %% result
# a dynamic graph, repeat playing
idx = 4
# each epoch costs "interval" millisecond
# my_plot2d_animate(states_epoch, idx=idx, interval=200, my_cmap='Greys')
my_plot2d_animate(states_epoch, idx=idx, interval=200)  # None=='viridis', **high contrast**
# my_plot2d_animate(states_epoch, idx=idx, interval=200, my_cmap='bg')  # brown green, project
# other choices: mpl.cm.bwr,
# %% for debugging purposes
time = 6
# time = -1
idx = 0
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
