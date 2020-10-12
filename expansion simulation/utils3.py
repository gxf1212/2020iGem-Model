import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

# we call the cell in CA "square" to distinguish from biological "cell"
global ELE
ELE = 6  # number of elements in a parallelogram


# %% define class cell
# class Cell(dtype='uint18'):
#     '''
#     # attributes:
#     b_BS, s_BS: the number of B. subtilis in two states
#     b_No, s_No: Nostoc sp.
#     N: nutrition content
#     EPS: EPS content
#     '''
#
#     def __init__(self):
#         self.b_BS=self.s_BS=self.b_No=self.s_No=0
#         self.N=0
#         self.EPS=0


# %% basic and common tool functions
def pairs(n: int):
    # return the indexes in the form of real number pairs
    # and this func is not considering the phases of new bacteria (put it later...)
    idx = range(n ** 2)  # random sequence
    row = np.array([int(i / n) for i in idx])
    col = np.array(idx) % n
    return [(row[i], col[i]) for i in range(len(idx))]


def rand_pairs(n: int):
    # random index pairs
    idx = np.random.permutation(n ** 2)  # random sequence
    row = np.array([int(i / n) for i in idx])
    col = np.array(idx) % n
    return [(row[i], col[i]) for i in range(len(idx))]


def find_neighbors(idx_pair, n, grid="rect_Moore"):
    # for rectangular, this is easy
    # we use Moore neighborhood/Von Neumann neighborhood
    x, y = idx_pair
    neighbors = []
    if grid == "rect_Moore":
        for i in [xi for xi in (x - 1, x, x + 1) if 0 <= xi < n]:  # exclude outer squares
            for j in [yi for yi in (y - 1, y, y + 1) if 0 <= yi < n]:
                neighbors.append((i, j))

    if grid == "rect_Neumann":
        for i in [xi for xi in (x - 1, x, x + 1) if 0 <= xi < n]:
            neighbors.append((i, y))
        for j in [yi for yi in (y - 1, y, y + 1) if 0 <= yi < n]:
            neighbors.append((x, j))
        neighbors.remove(idx_pair)  # idx_pair was added twice

    # for a hexagon grid, namely a rhombus
    # still using a matrix to store, but we're using just 3/4 of its area
    # the rhombus leans to the right, so we neglect the upper left and lower right
    if grid == "hexagonal":
        for i in [xi for xi in (x - 1, x, x + 1) if 0 <= xi < n]:  # exclude outer squares
            for j in [yi for yi in (y - 1, y, y + 1) if 0 <= yi < n]:
                neighbors.append((i, j))
        # exclude two extra parallelogram
        try:
            neighbors.remove((x - 1, y + 1))
        except:
            pass
        try:
            neighbors.remove((x + 1, y - 1))
        except:
            pass

    neighbors.remove(idx_pair)  # exclude itself

    return neighbors  # at most 8 neighbors around the square


def reverse_neighbor_value(idx_pair, nei, state_init):
    """
    it's wrong!
    diff.append(np.array([reverse_neighbor_value(idx_pair, nei, state_init[i]) + state_init[i][nei]
                    - 2 * state_init[i][idx_pair] for nei in neighbors], dtype='float64'))
    :param: state_init: of a particular element!
    :return: the number in reverse neighbor. if out of border, return 0
    just find the symmetry point
    """
    n = state_init.shape[-1]
    x, y = idx_pair
    xn, yn = nei
    xr = 2 * x - xn
    yr = 2 * y - yn
    if 0 <= int(xr) < n and 0 <= int(yr) < n:
        return state_init[xr][yr]
    else:
        return 0


# %% initialization methods
def init_classic_center(n: int, num=(100, 200, 100, 100), dtype='uint16'):
    # create initial distribution by putting some bacteria in the center
    # only "active" cell
    # B.S and Nostoc, respectively
    state0 = np.zeros((ELE, n, n), dtype=dtype)
    idx = [0, 2, 4, 5]
    for i in idx:
        state0[i, int(n / 2), int(n / 2)] = num[idx.index(i)]

    return state0


def init_free(n: int, location=(0, 0), num=(10, 4, 10, 4, 10, 10), dtype='uint16'):
    """
    create initial distribution by setting one parallelogram with certain state and number
    if you want to set >1 seeds, please add multiple objects this function returns with different params
    :param num: a tuple containing the numbers of each states
    """
    state0 = np.zeros((ELE, n, n), dtype=dtype)
    for i in range(len(num)):
        state0[i, location] = num[i]

    return state0


# %% numerical tool functions, e.g calculate gradients
def normalize(x):
    # normalize with mu and sigma
    mu = np.mean(x)
    sigma = np.mean(np.square(x - mu))
    if sigma == 0:
        sigma = 1
    result = (x - mu) / sigma

    return result


def calcu_diff(state_init, idx_pair, neighbors):
    """
    :param state_init: always the same in one epoch
    :param grid: for rectangular grid: Moore (8) or Neumann (4); or hexagonal (6)
    :return: the local 2nd derivative of BS, NO, EPS, nutrient
    spores are ignored because we ignore water flux
    nutrient != all molecules, but can be controlled by probability
    """
    var_idx = [0, 2, 4, 5]  # for BS, NO, EPS, nutrient
    if state_init.ndim == 2:  # a single microorganism
        var_idx = [0]
        state_init = [state_init]  # add a dimension
    diff = []
    for i in var_idx:  # equal to original [0,2,4,5]
        # get differences in neighbors' order. convert to float in case of negative numbers and following calculation
        diff.append(np.array([state_init[i][idx_pair] - state_init[i][nei] for nei in neighbors], dtype='float64'))
    return np.array(diff)  # a row vector for each variable, shape = (len(var_idx), len(neighbors))


# %% migration
def calcu_num_migration(diff, weight, probs_migration, dtype):
    """
    calculate number of cell migration. cps is calculated in function "update_number_mig"
    TODO: if we decide to ignore water flux, those commented statements will be deleted in this function and update
    :param diff: 2nd derivative. 0~3 for BS, NO, EPS, nutrient. see calcu_diff
    :param weight: 2 tuples (each len=3). weight for BS, No, EPS and nutrient.
    set according to the degree of impact. EPS and nutrient are not spreading by itself
    :param probs_migration: list. len=4. [2,3] is now unused
    migration probabilities, actually diffusion coefficient of BS, No, EPS and nutrient
    :param dtype: we must control dtype here
    :return: diffusion coefficient * 2nd derivatives * normalized weights
    due to 2nd derivative assumption, we are not multiplying states[(i, 2), x, y].reshape(2, 1)
    we add because bacteria both move themselves and take away by osmotic pressure & swelling at the same time
    weights contain info on both the importance and unit unification between EPS and bacteria
    after adding all mechanisms, we just use the positive term, because num_migration in opposite
    direction are opposite numbers
   """
    num_migration = np.zeros(shape=diff.shape, dtype=dtype)
    for i in [0, 1]:  # BS, No, a row for all neighbors
        wei = np.array([weight[i]]).T
        temp = probs_migration[i] * diff[[i, 2, 3], :] * wei / np.sum(wei)
        num_migration[i, :] = np.sum(temp, axis=0).round()
    # for i in [2, 3]:
    #     num_migration[i, :] = probs_migration[i] * diff[i, :]
    num_migration = np.maximum(num_migration, 0)  # discard negative terms

    return num_migration


def update_number_mig(state_update, idx_pair, neighbors, num_migration, ratio_cps):
    """
    update population numbers (only once for a certain parallelogram)
    CPS will be taken away by microorganism but nutrient won't
    :param ratio_cps: how much units of CPS a BS/No can take away (or, how much a cell contains)
    """
    # idx = [0,2,4,5]
    idx = [0, 2]
    # update the neighbors
    for nei in neighbors:
        # cell migration
        for i in idx:
            # in case num_mig > num_center. directly change the values in this array
            num_migration[idx.index(i), neighbors.index(nei)] = np.minimum(state_update[i][idx_pair], num_migration[idx.index(i), neighbors.index(nei)])
            state_update[i][nei] += num_migration[idx.index(i), neighbors.index(nei)]
            state_update[i][idx_pair] -= num_migration[idx.index(i), neighbors.index(nei)]

        # calculate cps migration, using modified num_migration
        for i in [0, 1]:
            state_update[4][nei] += num_migration[i, neighbors.index(nei)] * ratio_cps[i]

    return state_update


def update_mig(idx_pairs, n, state_update, state_init, grid,
               probs_migration, weight, ratio_cps,
               dtype='uint16'):
    """
    rectangular/hexagonal share this function
    :param idx_pairs: sequential pairs
    :param weight: to determine the effect of population, EPS, nutrient
    :param probs_migration: diffusion coefficient. as a intrinsic property of microorganisms and substances
    TODO: to be unified into this CA unit
    :var: idx_pair: the center
    :return: new state (t+1/2 in the paper)
    based on diffusion equation, where particles move according to gradient
    we just consider the difference between the center and its neighbors to determine the 2nd derivative
    TODO: may be changed, using random variables. now we round off the numbers
    before doing so, it should all be isotropic
    """
    for idx_pair in idx_pairs:
        neighbors = find_neighbors(idx_pair, n, grid)  # find neighbors
        # if state_init[0][idx_pair] != 0:  # for debugging, find where there are bacteria
        #     print("not empty")
        diff = calcu_diff(state_init, idx_pair, neighbors)
        num_migration = calcu_num_migration(diff, weight, probs_migration, dtype)
        state_update = update_number_mig(state_update, idx_pair, neighbors, num_migration, ratio_cps)

    return state_update


#%% consolidation
def going_out(idx_pair, neighbors, num, state_init, state_update):
    """
    redundant cells go out, the mechanism of cell "push"!
    we assume that this push follows the gradient of cell density, not simply outwards
    :param idx_pair: center
    :param num: to go out
    :param state_init, state_update are both states of a single microorganism
    :return: new state (set state[idx_pair] in the parent function)
    """
    diff = calcu_diff(state_init, idx_pair, neighbors).squeeze()  # generally they will all be positive
    diff = np.maximum(diff, 0)  # ignore negative migration. shape = (1, len(neighbors))
    num_out = (diff / np.sum(diff) * num).round()  # distribute num according to proportion of diff
    # print(diff)
    print(num_out)
    for nei in neighbors:
        state_update[nei] = state_update[nei] + num_out[neighbors.index(nei)]  # add neighbors

    return state_update


def update_con_grow(n, grid, idx_pairs_rand, r, K, state_update, ratio=2):
    """
    :param state_update: only input number of BS (1) and No (2)
    :param r: a tuple of (r1, r2)
    :param K: a tuple of (K1, K2)
    :param ratio: a certain ratio of K to limit growth. 2 (default) has no effect
    :return new state of BS and No
    if num < K, just grow; or the redundant cells are pushed out ("outwards?")
    """
    # grow
    for i in [0, 1]:
        grow = np.ones(shape=state_update[i].shape) * (1 + r[i])  # (n,n), times to multiply
        for idx_pair in idx_pairs_rand:
            if state_update[i][idx_pair] > ratio * K[i]:  # number > a certain ratio of K
                grow[idx_pair] = 1  # cell won't grow
        state_update[i] = state_update[i] * grow

    state_init = state_update  # copy, for going through idx pairs
    # expansion due to overgrowth and redundant cells
    while True:  # iterate until all parallelograms' values < K[i]
        for i in [0, 1]:
            big = [idx_pair for idx_pair in idx_pairs_rand if state_init[i][idx_pair] > K[i]]  # find the bigger idxes
            for idx_pair in big:  # if big, some goes out
                neighbors = find_neighbors(idx_pair, n, grid)  # find neighbors
                state_update[i] = going_out(idx_pair, neighbors, state_update[i][idx_pair] - K[i],
                                            state_init[i], state_update[i])
                state_update[i][idx_pair] = K[i]

            if np.max(state_update[0]) <= K[0] and np.max(state_update[1]) <= K[1]:  # judge the two simultaneously
                # may return going_out to calculate CPS migration
                return state_update[0], state_update[1]
            # print("still iterating")  # for debugging


def update_con_spores(state_update, m1, m2):
    # we ignore natural decay
    # spores are only caused by sufficient nutrient

    return state_update


def update_con_nutrient(state_update):
    return state_update


def update_con(idx_pairs_rand, n, state_update, grid, dtype,
               params_BS, params_No):
    """
    1 for BS, 2 for No
    :param r1: birth rate. alpha
    :param l1: turn into spores. relatively low
    :param y1: nutrition consumption rate. gamma in the paper. we are not separating Q and respiratory rate
    :return: state_update, new state (b^t+1)
    we assume that the state has something to do with update order, because the process of pushing.
    this order should be random
    """

    r1, l1, y1, K1 = params_BS
    r2, l2, y2, K2 = params_No

    state_init = state_update.copy()
    # if state_init[0][idx_pair] != 0:  # for debugging, find where there are bacteria
    #     print("not empty")
    # grow; redundant cells going out. only use states of BS and No
    state_update[0], state_update[2] = update_con_grow(n, grid, idx_pairs_rand,
                                                       (r1, r2), (K1, K2), state_update[[0, 2]], 2)
    # after growth, some turn into spores
    # state_update = update_con_spores()
    # consume / produce nutrients (carbon), using the number of new cells
    # state_update = update_con_nutrient()

    return state_update


# %% main function
def stimulation_v3(n, grid='rect_Moore', state_init=None, epoch=10,
                   see_phase=False, dtype='uint16',
                   probs_migration=(0.1, 0.1, 0.1, 0.1), weight=((1, 1, 1), (1, 1, 1), 1, 1),
                   ratio_cps=(1, 1),
                   params_BS=(0.5, 0.05, 1, 0.01, 600, 20), params_No=(0.5, 0.05, 1, 0.01, 600, 20)):
    """
    we call the overall cycle/period "epoch"
    we call the two phases (migration; consolidation) "phase"
    each phase involves multiple "updates", going through all points and values
    :param grid: for rectangular grid: Moore (8) or Neumann (4); or hexagonal (6)
    :param n:
    :param grid:
    :param state_init:
    :param epoch:
    :param see_phase:
    :param dtype:
    :param probs_migration:
    :param weight:
    :param ratio_cps:
    :param params_BS:
    :param params_No:
    :return:
    """

    if state_init is None:
        print("Error!")
        return

    # preparation
    states_epoch = np.zeros(shape=(epoch + 1, ELE, n, n), dtype=dtype)
    state_update = states_epoch[0] = state_init  # create a 3 dimensional array and initialize
    idx_pairs = pairs(n)  # generate sequential pairs of location coordination
    idx_pairs_rand = rand_pairs(n)
    if see_phase:  # only if we want to observe the change after each phase, we add this, or it costs time
        states_phase = np.zeros(shape=(epoch * 2 + 1, ELE, n, n), dtype=dtype)
    else:
        states_phase = None

    # simulation
    for time in range(epoch):
        """
        state_init: keep the state at the start of this phase, and use its values
        state_update: state after each update, continuously changing
        copy means a no relationship with the original one
        """
        # migration phase

        state_update = update_mig(idx_pairs, n,
                                  state_update, state_init=state_update.copy(), grid=grid,
                                  probs_migration=probs_migration, weight=weight, ratio_cps=ratio_cps,
                                  dtype=dtype)

        if see_phase:
            states_phase[time * 2 + 1] = state_update.copy()

        # consolidation phase
        state_update = update_con(idx_pairs_rand, n,
                                  state_update, grid=grid,
                                  dtype=dtype, params_BS=params_BS, params_No=params_No)

        if see_phase:
            states_phase[time * 2 + 2] = state_update.copy()

        states_epoch[time + 1] = state_update
        print("The iteration epoch "+str(time)+" has finished normally")
    return states_epoch, states_phase


# %% visualization
# see https://matplotlib.org/tutorials/colors/colormaps.html for more colormaps
def my_plot2d_animate(ca, idx=0, title='evolve', interval=50, my_cmap='Greys'):
    cmap = plt.get_cmap(my_cmap)  # color
    fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(1, 1, 1)
    im = plt.imshow(ca[0][idx], animated=True, cmap=cmap)  # set the value range here!, vmin=0, vmax=3
    plt.colorbar(ticks=get_ticks(np.max(ca)))  # set special values in the bar
    i = {'index': 0}  # because it is used in a newly-defined function

    def updatefig(*args):
        i['index'] += 1
        if i['index'] == len(ca):
            i['index'] = 0
        title = 'evolve at epoch ' + str(i['index'])  # manually set 'epoch' or 'stage'
        ax.set_title(title)
        im.set_array(ca[i['index']][idx])
        return im,

    ani = animation.FuncAnimation(fig, updatefig, interval=interval, blit=True)
    # ani.save(filename="evolve.gif",writer='pillow')  # save as a gif file, change writer thus no error
    plt.show()


def my_plot2d(ca, timestep=None, title='', my_cmap='Greys', idx=0):
    cmap = plt.get_cmap(my_cmap)
    fig = plt.figure(figsize=(9.6, 7.2))
    plt.title(title)
    if timestep is None:  # maybe the input is a single-time state
        data = ca[idx]
    else:
        data = ca[timestep, idx]
    plt.imshow(data, interpolation='none', cmap=cmap)
    plt.colorbar(ticks=get_ticks(np.max(data)))


def get_ticks(max, sep=1000):
    return list(sep * np.arange(int(max / sep) + 1))+[max]
