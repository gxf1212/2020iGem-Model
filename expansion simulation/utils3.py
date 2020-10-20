import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

# we call the cell in CA "square" to distinguish from biological "cell"
global ELE
ELE = 6  # number of elements in a parallelogram
dict_cell = {0: "BS", 1: "Nostoc"}  # for outputting information

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
def init_classic_center(n: int, num=(500, 500, 200, 1000), dtype='uint16'):
    """
    create initial distribution by putting some bacteria in the center
    :param num: B.S and Nostoc, respectively. only "active" cell.
    note: initial EPS and nutrient is put everywhere with the value given
    guarantee initial value is enough by adding cps at each seed
    :return: initial state
    """
    state0 = np.zeros((ELE, n, n), dtype=dtype)
    for i in [0, 2]:
        state0[i, int(n / 2), int(n / 2)] = num[int(i/2)]
    for i in [4, 5]:
        state0[i] = np.ones(shape=(n, n), dtype=dtype) * num[i-2]
    state0[4, int(n / 2), int(n / 2)] += num[0] + num[1]

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


def init_multi_pos(n:int, idx_pairs, cell_values, eps, nut, dtype='uint16'):
    """
    :param idx_pairs: tuple (2). all positions you want to place bacteria
    :param cell_values: tuple (2). according to idx_pairs, BS and No values
    :param eps, nut: for everywhere
    :return: initial state
    """
    state0 = np.zeros((ELE, n, n), dtype=dtype)
    # bacteria
    for i in range(len(idx_pairs)):
        for j in [0, 1]:
            state0[j][idx_pairs[i]] = cell_values[i][j]
            state0[4][idx_pairs[i]] += cell_values[i][j]
    # eps and nutrient
    state0 += init_classic_center(n, (0, 0, eps, nut), dtype)

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
    due to gradient assumption, we are not multiplying states[(i, 2), x, y].reshape(2, 1)
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


def update_cps(state_start, state_end, ratio_cps):
    """
    called by pure migration function: update_number_mig, update_con_grow/push
    we just calculate difference between cell numbers, and * cps ratio
    :param state_start: start
    :param state_end: end
    :param ratio_cps: 0 for BS, 1 for No
    :return: new state(end, update), only updated eps
    """
    state_update = state_end
    for i in [0, 1]:
        diff = state_end[i] - state_start[i]  # >0 means more cps
        state_update[4] += diff * ratio_cps[i]

    if np.min(state_update[4]) < 0:  # no enough eps to take away
        return state_update, True
    else:
        return state_update, False


def update_number_mig(state_update, idx_pair, neighbors, num_migration, ratio_cps):
    """
    update population numbers (only once for a certain parallelogram)
    CPS will be taken away by microorganism but nutrient won't
    there might be much more EPS units, but a unit is small. we can guarantee there are enough EPS for bacteria to take with
    :param ratio_cps: how much units of CPS a BS/No can take away (or, how much a cell contains)
    """
    # idx = [0,2,4,5]
    idx = [0, 2]
    # cell migration, driven by osmotic pressure
    for nei in neighbors:
        for i in idx:
            # in case num_mig > num_center. directly change the values in this array
            num_migration[idx.index(i), neighbors.index(nei)] = np.minimum(state_update[i][idx_pair], num_migration[idx.index(i), neighbors.index(nei)])
            state_update[i][nei] += num_migration[idx.index(i), neighbors.index(nei)]
            state_update[i][idx_pair] -= num_migration[idx.index(i), neighbors.index(nei)]

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
    - based on diffusion equation, where particles move according to gradient
    we just consider the difference between the center and its neighbors to determine the 2nd derivative
    TODO: may be changed, using random variables. now we round off the numbers
    before doing so, it should all be isotropic
    """
    state_start = state_update.copy()  # for updating cps
    for idx_pair in idx_pairs:
        neighbors = find_neighbors(idx_pair, n, grid)  # find neighbors
        # if state_init[0][idx_pair] != 0:  # for debugging, find where there are bacteria
        #     print("not empty")
        diff = calcu_diff(state_init, idx_pair, neighbors)
        num_migration = calcu_num_migration(diff, weight, probs_migration, dtype)
        state_update = update_number_mig(state_update, idx_pair, neighbors, num_migration, ratio_cps)

    # cps migration, taken by cells, using modified num_migration
    state_update, exit_flag = update_cps(state_start, state_update.copy(), ratio_cps)
    if exit_flag:
        print("CPS is not enough when updating migration.\nProgram will exit.")

    return state_update, exit_flag


#%% consolidation
def going_out(idx_pair, neighbors, num, state_update):
    """
    redundant cells go out, the mechanism of cell "push"!
    we assume that this push follows the gradient of cell density, not simply outwards
    but we just round them off
    :param idx_pair: center
    :param num: to go out
    :param state_update is states of a single microorganism
    :return: new state (set state[idx_pair] in the parent function)
    """
    diff = calcu_diff(state_update, idx_pair, neighbors).squeeze()  # generally they will all be positive
    diff = np.maximum(diff, 0)  # ignore negative migration. shape = (1, len(neighbors))
    if np.max(diff) == 0:  # if diff are all zero
        return state_update
    num_out = (diff / np.sum(diff) * num).round()  # distribute num according to proportion of diff
    for nei in neighbors:
        try:
            state_update[nei] = state_update[nei] + num_out[neighbors.index(nei)]  # add neighbors
        except Exception as e:
            print(num_out[neighbors.index(nei)])

    return state_update


def update_con_grow(idx_pairs_rand, r, K, state_update, ratio=2):
    """
    growth. the numbers increase
    :param state_update: only input number of BS (1) and No (2)
    :param r: a tuple of (r1, r2)
    :param K: a tuple of (K1, K2)
    :param ratio: a certain ratio of K to limit growth. 2 (default) has no effect
    :return new state of BS and No
    """
    for i in [0, 1]:
        grow = np.ones(shape=state_update[i].shape) * (1 + r[i])  # (n,n), times to multiply
        for idx_pair in idx_pairs_rand:
            if state_update[i][idx_pair] > ratio * K[i]:  # number > a certain ratio of K
                grow[idx_pair] = 1  # cell won't grow
        state_update[i] = state_update[i] * grow
        if np.sum(state_update[i]) / state_update.shape[-1] / state_update.shape[-2] > K[i]:
            print("The total number of " + dict_cell[i] + " cells has exceeded the total carrying capacity. "
                  "\nNo going-out will achieve a balance. Please increase n. Program will exit.")
            # return exit = True to go back to function "simulation" and output
            return state_update[0], state_update[1], True
        else:
            return state_update[0], state_update[1], False


def update_con_push(n, grid, idx_pairs_rand, K, state_update):
    """
    expansion due to overgrowth and redundant cells. if not so, just return back
    :param state_update: only input number of BS (1) and No (2)
    :return new state of BS and No
    if num < K, just grow; or the redundant cells are pushed out ("outwards?")
    - we assume that cells go out in a random order, because the process of pushing.
    so we just use state_update to calculate diff...just add state_init if you want to update simultaneously
    - we assume that cells are pushed based merely on cell density gradient
    TODO: may add "action at a distance" (超距作用), if some consecutive cells exceed K, they just push the outest one
    to replace max_iter
    """
    max_iter = 100  # if the cells are crowded, they need a few iterations to achieve <= K
    for t in range(max_iter):  # try to iterate until all parallelograms' values < K[i]
        for i in [0, 1]:
            big = [idx_pair for idx_pair in idx_pairs_rand if state_update[i][idx_pair] > K[i]]  # find the bigger idxes
            for idx_pair in big:  # if big, some goes out
                neighbors = find_neighbors(idx_pair, n, grid)  # find neighbors
                state_update[i] = going_out(idx_pair, neighbors, state_update[i][idx_pair] - K[i],
                                            state_update[i])
                state_update[i][idx_pair] = K[i]

            if np.max(state_update[0]) <= K[0] and np.max(state_update[1]) <= K[1]:  # judge the two simultaneously
                # may return going_out to calculate CPS migration
                return state_update[0], state_update[1], False
            # print("iteration" + str(t))  # for debugging

    print("After " + str(max_iter) + " iterations' pushing, the cells are still too dense "
                                     "to push each other and reach K"
                                     "\nProgram will exit.")

    return state_update[0], state_update[1], True


def update_con_nutrient(state_update, l, y, p, ratio_BN, dtype):
    """
    we ignore natural decay
    only sufficient nutrient cause spores
    :param state_update:
    :param l: a tuple of (l1, l2)
    :param y: a tuple of (y1, y2)
    :param p: a tuple of (p1, p2)
    :param ratio_BN: determining how much N/P BS can provide for No
    :return: new state of all
    a unit nutrient means the least amount for a BS to consume in a single timestep
    we assume nutrient is the limiting factor of BS, N and P (number of BS) for No
    """
    # determine nutritional condition
    # if nutrient is not enough to support bacteria, those redundant ones turn into spores
    # those positions are False, thus not participating the following calculation
    num_spo1 = (state_update[5] < state_update[0]) * (state_update[0] - state_update[5])
    state_update[0] -= num_spo1
    state_update[1] += num_spo1
    # if nutrient is much more, just recover
    num_rec1 = (state_update[5] > state_update[0] + state_update[1]) * np.array(state_update[1] * l[0], dtype)
    state_update[1] -= num_rec1
    state_update[0] += num_rec1
    # the same for No, but rely on a ratio
    # if BS is not enough to provide N/P, too low ratio or zero
    # state_update[0] == 0 or
    num_spo2 = (state_update[0] < state_update[2] * ratio_BN)\
               * np.array(state_update[2] - state_update[0] / ratio_BN, dtype)
    state_update[2] -= num_spo2
    state_update[3] += num_spo2
    # if BS is enough
    # or state_update[0] == 0
    num_rec2 = (state_update[0] > (state_update[2] + state_update[3]) * ratio_BN)\
               * np.array(state_update[3] * l[1], dtype)
    state_update[3] -= num_rec2
    state_update[2] += num_rec2

    # update nutrient. minus. y[No] < 0
    for i in [0, 1]:
        state_update[5] -= np.array(state_update[i] * y[i], dtype)

    # update rps
    for i in [0, 1]:
        state_update[4] += np.array(state_update[i*2] * p[i], dtype)  # s[2]*p[1]

    return state_update


def update_con(idx_pairs_rand, n, state_update, grid, dtype,
               params_BS, params_No, ratio_cps, ratio_BN=1):
    """
    1 for BS, 2 for No
    :param r: birth rate. alpha
    :param l: rate of change into spores when nutrients are enough. maybe include death due to slow recovery
    :param y: nutrition consumption rate. gamma in the paper. we are not separating Q and respiratory rate
    :param p: rps production rate. not too big (depending on ratio of cps and rps)
    :return: state_update, new state (b^t+1)
    """
    r1, l1, y1, p1, K1 = params_BS
    r2, l2, y2, p2, K2 = params_No

    # if state_init[0][idx_pair] != 0:  # for debugging, find where there are bacteria
    #     print("not empty")

    # consume / produce nutrients (carbon) first. some turn into spores after migration and nutrient consumption
    state_update = update_con_nutrient(state_update, (l1, l2), (y1, y2), (p1, p2), ratio_BN, dtype)

    state_start = state_update.copy()  # for calculating cps
    # grow; redundant cells going out. only use states of BS and No
    state_update[0], state_update[2], exit_flag = update_con_grow(idx_pairs_rand, (r1, r2),
                                                                  (K1, K2), state_update[[0, 2]], 2)
    if exit_flag:  # if exit == True, stop the simulation immediately, pass the flag to "simulation"
        state_update, exit_flag2 = update_cps(state_start, state_update, ratio_cps)
        if exit_flag2:
            print("CPS is not enough when updating growing.\nProgram will exit.")
        return state_update, exit_flag

    state_update[0], state_update[2], exit_flag = update_con_push(n, grid, idx_pairs_rand,
                                                                  (K1, K2), state_update[[0, 2]])
    state_update, exit_flag2 = update_cps(state_start, state_update.copy(), ratio_cps)
    if exit_flag2:
        print("CPS is not enough when updating pushing.\nProgram will exit.")
    # no matter exit_flag is, it will return and update cps. but only if flags are all false, return false

    return state_update, exit_flag or exit_flag2


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
    :param n: size of grid
    :param state_init: we set one / multiple seeds
    :param epoch:
    :param see_phase: if we return states after the migration phase (i+1/2 epochs)
    :param dtype: 'int16'
    :param probs_migration:
    :param weight:
    :param ratio_cps:
    :param params_BS:
    :param params_No:
    :return: states along the time
    notes:
    -
    """

    if state_init is None:
        print("Error!")
        return

    # preparation
    states_epoch = [state_init]
    state_update = state_init  # create a 3 dimensional array and initialize
    idx_pairs = pairs(n)  # generate sequential pairs of location coordination
    idx_pairs_rand = rand_pairs(n)
    states_phase = []

    # simulation
    for time in range(epoch):
        """
        state_init: keep the state at the start of this phase, and use its values
        state_update: state after each update, continuously changing
        copy means a no relationship with the original one
        """
        # migration phase
        state_update, exit_flag = update_mig(idx_pairs, n,
                                             state_update, state_init=state_update.copy(), grid=grid,
                                             probs_migration=probs_migration, weight=weight, ratio_cps=ratio_cps,
                                             dtype=dtype)
        if exit_flag:  # stop the simulation immediately
            return np.array(states_epoch, dtype=dtype), np.array(states_phase, dtype=dtype)

        if see_phase:  # only if we want to observe the change after each phase, we add this, or it will cost time
            states_phase.append(state_update.copy())

        # consolidation phase
        state_update, exit_flag = update_con(idx_pairs_rand, n,
                                             state_update, grid=grid, dtype=dtype,
                                             params_BS=params_BS, params_No=params_No, ratio_cps=ratio_cps)
        if exit_flag:  # stop the simulation immediately
            return np.array(states_epoch, dtype=dtype), np.array(states_phase, dtype=dtype)

        if see_phase:
            states_phase.append(state_update.copy())

        states_epoch.append(state_update)
        print("The iteration epoch "+str(time)+" has finished normally")

    return np.array(states_epoch, dtype=dtype), np.array(states_phase, dtype=dtype)


# %% visualization
# see https://matplotlib.org/tutorials/colors/colormaps.html for more colormaps
def my_plot2d_animate(ca, idx=0, title='evolve', interval=50, my_cmap='Greys'):
    cmap = plt.get_cmap(my_cmap)  # color
    fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(1, 1, 1)
    im = plt.imshow(ca[0][idx], animated=True, cmap=cmap)  # set the value range here!, vmin=0, vmax=3
    plt.colorbar(ticks=get_ticks(np.max(ca[:, idx])))  # set special values in the bar
    i = {'index': 0}  # because it is used in a newly-defined function

    def updatefig(*args):
        i['index'] += 1
        if i['index'] == len(ca):
            i['index'] = 0
        title = 'evolve at epoch ' + str(i['index']) + ' for index ' + str(idx)  # manually set 'epoch' or 'stage'
        ax.set_title(title)
        im.set_array(ca[i['index']][idx])
        return im,

    ani = animation.FuncAnimation(fig, updatefig, interval=interval, blit=True)
    # ani.save(filename="evolve.gif",writer='pillow')  # save as a gif file, change writer thus no error
    plt.show()


def my_plot2d(ca, timestep=None, my_cmap='Greys', idx=0):
    cmap = plt.get_cmap(my_cmap)
    fig = plt.figure(figsize=(9.6, 7.2))
    if timestep is None:  # maybe the input is a single-time state
        data = ca[idx]
        title = 'evolve for index ' + str(idx)
    else:
        data = ca[timestep][idx]
        title = 'evolve at epoch ' + str(timestep) + ' for index ' + str(idx)
    plt.title(title)
    plt.imshow(data, interpolation='none', cmap=cmap)
    plt.colorbar(ticks=get_ticks(np.max(data)))


def get_ticks(max, sep=None):
    if sep is None:
        if max < 600:
            sep = 100
        else:
            if max < 1200:
                sep = 200
            else:
                sep = 1000
    ticks = list(sep * np.arange(int(max / sep) + 1))
    if max > ticks[-1] + 0.1 * sep:  # to avoid the situation that max is a little bigger than the biggest tick
        ticks = ticks + [max]
    return ticks  # [0] + ?
