# try simple CA
# we use rectangular grid
# dtype='uint8'

from utils12 import *
import cellpylib as cpl

# %% parameters
# fundamental parameters
'''type of states
0 for empty
1 for migration, move to a neighborhood "cell"
2 for proliferation, put one cell in randomly direction
3 for quiescence, no action
# 4 for death
'''
n = 200  # size
Tp1 = 1  # cell 2 -> 2, 3 is so big that program exits early
Tp2 = 1  # cell 2 -> 1
Tm = 2  # migration

# state transfer matrix, trans_1 stands for probability of changing from state 1 to 1,2,3
trans_1 = (2, 7, 1)
trans_2 = (0, 1, 0)
trans_3 = (0, 0, 1)

epoch = 20  # number of period

# %% easy stimulation
state_init = cpl.init_simple2d(n, n, val=2, dtype='uint8')
states = stimulation_v1(n, state_init, epoch=epoch)

# %% result
# cpl.plot2d(states, timestep=1)  # visualize single timestep
# cpl.plot2d_animate(states)  # a dynamic graph
# np.mean(states[0]==states[5])

# %% 3-stage rule with single BS and no nutrition
# the init_simple2d added a dimension to a state, which fits the original evolution function
# put a state 2 cell in the center, but we can write our own
# each time I run the simulation, state_init changed ????! have to put it here
state_init = cpl.init_simple2d(n, n, val=2, dtype='uint8').squeeze()  #
states_epoch, states_stage = stimulation_v2(n=n, state_init=state_init,
                                                   Tp1=Tp1, Tp2=Tp2, Tm=Tm, epoch=epoch,
                                                   see_stage=True, trans_1=trans_1)
# it is still hard to perform a large-scale simulation
'''
下次要讨论扩散规则！！！交换想法，梳理我的，明确目标。确定讨论内容，包括生长参数的
他们本身扩散欲望不强（那些多糖会帮助粘在沙粒上，念珠藻还用藻丝纤维缠绕固定土壤，当然有一部分是分泌的），只是因为长满才扩散。
如果结合实际，我们应该找到使二者能够扩散的因素，比如，就因为长满才扩散。现在的程序和框架可用，但这算是和实际结合的一个点
（和多糖平行，大家考虑）
'''
# %% result
# a dynamic graph, repeat playing
my_plot2d_animate(states_epoch, interval=600)  # each epoch costs 800 millisecond
# visualize a single timestep
# my_plot2d(states_epoch, timestep=1)


# %%
# for debugging purposes
my_plot2d_animate(states_stage)  # dynamic graph by stage
# my_plot2d(states_stage, timestep=2)

# my_plot2d(state_update)
# my_plot2d(state_init)
# my_plot2d(state_init_round)
# my_plot2d(states_epoch, timestep=-1)
# np.where(states_epoch[-1])==3
# np.where(state_init_round==2)  # np.where() find some certain element in a state, returning the index
