# try CA version 2, Modelling the morphology of migrating bacterial colonies
# we use rectangular grid or hexagonal grid
# dtype='uint16', 0~65535

from utils3 import *

global ELE
dtype='int16'


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
epoch = 50  # number of period. more causes NaN??
probs_migration = [0.02, 0.02, 0.02, 0.02]  # diffusion coefficient
weight = [[1, 1, 1], [1, 1, 1]]  # cell, EPS and nutrient on migration
params_BS = (0.1, 0.03, 0.1, 0.5, 1000)  # growth, to spores, nutrient consumption, rps production, carrying capacity
params_No = (0.1, 0.03, -0.1, 0.5, 1000)

# the init_simple2d added a dimension to a state, which fits the original evolution function
# put a state 2 cell in the center, but we can write our own
# each time I run the simulation, state_init changed ????! have to put it here
state_init = init_classic_center(n, num=(200, 200, 200, 1000), dtype=dtype)
states_epoch, states_phase = stimulation_v3(n=n, state_init=state_init, grid='rect_Moore',
                                            # rect_Moore, rect_Neumann, hexagonal
                                            epoch=epoch, see_phase=False,
                                            probs_migration=probs_migration, weight=weight,
                                            params_BS=params_BS, params_No=params_No,
                                            dtype=dtype)
# it is still hard to perform a large-scale simulation

# %% result
# a dynamic graph, repeat playing
idx = 3
# my_plot2d_animate(states_epoch, idx=idx, interval=600)  # each epoch costs "interval" millisecond
my_plot2d_animate(states_epoch, idx=idx, interval=400)


# %% for debugging purposes

# time = 20
time = -1
idx = 4
my_plot2d(states_epoch, timestep=time, idx=idx)
# states_epoch[time]
states_epoch[time][idx]
# my_plot2d(state_update, idx=0)
# my_plot2d_animate(states_phase)  # dynamic graph by phase
# my_plot2d(states_phase, timestep=2)
# my_plot2d(state_init, idx=idx)
