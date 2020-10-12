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
n = 10  # size
epoch = 18  # number of period. more causes NaN??
probs_migration = [0.02, 0.02, 0.02, 0.02]  # diffusion coefficient
weight = [[1, 1, 1], [1, 1, 1]]  # cell, EPS and nutrient on migration
params_BS = (0.3, 0.03, 0.01, 1000)  # growth, to spores, nutrient consumption, carrying capacity
params_No = (0.3, 0.03, 0.01, 1000)

# the init_simple2d added a dimension to a state, which fits the original evolution function
# put a state 2 cell in the center, but we can write our own
# each time I run the simulation, state_init changed ????! have to put it here
state_init = init_classic_center(n, num=(1000, 1000, 1000, 1000), dtype=dtype)
states_epoch, states_phase = stimulation_v3(n=n, state_init=state_init, grid='rect_Moore',
                                            epoch=epoch, see_phase=False,
                                            probs_migration=probs_migration, weight=weight,
                                            params_BS=params_BS, params_No=params_No,
                                            dtype=dtype)
# it is still hard to perform a large-scale simulation
'''
交换想法，明确目标
他们本身扩散欲望不强（那些多糖会帮助粘在沙粒上，念珠藻还用藻丝纤维缠绕固定土壤，当然有一部分是分泌的），只是因为长满才扩散。
如果结合实际，我们应该找到使二者能够扩散的因素，比如，就因为长满才扩散。现在的程序和框架可用，但这算是和实际结合的一个点
（和多糖平行，大家考虑）
'''
# %% result
# a dynamic graph, repeat playing
idx = 2
# my_plot2d_animate(states_epoch, idx=idx, interval=600)  # each epoch costs "interval" millisecond
my_plot2d_animate(states_epoch, idx=idx, interval=600)


# %%
# for debugging purposes

time = 16
idx = 2
my_plot2d(states_epoch, timestep=time, idx=idx)
# states_epoch[time]
states_epoch[time][idx]
# my_plot2d(state_update, idx=0)
# my_plot2d_animate(states_phase)  # dynamic graph by phase
# my_plot2d(states_phase, timestep=2)

