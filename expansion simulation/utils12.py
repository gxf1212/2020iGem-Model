import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


# we call the cell in CA "square" to distinguish from biological "cell"


# %% basic and common tool functions
def rand_pair(n: int):
    # for rectangular, this is easy
    # return the indexes in the form of real number pairs
    # and this func is not considering the phases of new bacteria (put it later...)
    idx = np.random.permutation(n ** 2)  # random sequence
    row = np.array([int(i / n) for i in idx])
    col = np.array(idx) % n
    return [(row[i], col[i]) for i in range(len(idx))]


def find_neighbors_rect(idx_pair, n):
    # for rectangular, this is easy
    x, y = idx_pair
    neighbors = []
    for i in [xi for xi in (x - 1, x, x + 1) if 0 <= xi < n]:  # exclude outer squares
        for j in [yi for yi in (y - 1, y, y + 1) if 0 <= yi < n]:
            neighbors.append((i, j))
    neighbors.remove(idx_pair)  # exclude itself
    return neighbors  # at most 8 neighbors around the square


# %% initialization methods
def init_random(n: int, points=(0, 4, 0), dtype='uint8'):
    """
    create initial distribution by randomly setting some of the squares with certain state
    :param points: a tuple containing the numbers of each states
    """
    state0 = np.zeros((n, n), dtype=dtype)
    s1, s2, s3 = points
    state0[int(n / 2), int(n / 2)] = points

    return state0


# %% update and stimulation functions for previous versions
def update_easy_single(state, idx_pair, n: int):
    # update single point's state
    state_new = state
    neighbors = find_neighbors_rect(idx_pair, n)
    # find target for migration or proliferation. there must be more than 1. pick the 1st, which returns x,y
    target = np.random.permutation(len(neighbors))[0]
    update = False  # which can be substituted by using the state at the start of this stage
    if state[idx_pair] == 1:  # migration state
        state_new[idx_pair] = 0
        x, y = neighbors[target]
        state_new[x, y] = 1
        update = True

    return state_new, update


def stimulation_v1(n: int, state_init, epoch=20):
    states = np.zeros(shape=(epoch + 1, n, n), dtype='uint8')
    states[0] = state_init
    for t in range(epoch):
        idx_pairs = rand_pair(n)
        # actually it might be better to keep the order during each update
        for idx_pair in idx_pairs:
            states[t + 1], update = update_easy_single(states[t], idx_pair, n)  # not considering consecutive change!
            if update == True:
                break

    return states


# the above two is for test


def is_outwards(idx_pair1, idx_pair2, n):  # judge if point2 is on the outside of point1
    # now do nothing
    return True


def update_Tp1(state, idx_pair, n, option=8):  # 2 to 2
    state_new = state  # also, keep the state the start of this stage, for identifying 0/1/2
    neighbors = find_neighbors_rect(idx_pair, n)
    # default: use all neighbors. option can be adjusted
    # random indexes should not exceed len(neighbors)
    idx = np.random.permutation(len(neighbors))[:np.minimum(option, len(neighbors))]
    for i in idx:  # neighbors with selected indexes
        neighbor = neighbors[i]  # chosen neighbors
        if state[neighbor] == 0:  # if not empty (of its species), else no operation
            state_new[neighbor] = 2

    return state_new


def update_Tp2(state, idx_pair, n, option=8):  # 2 to 1
    state_new = state
    neighbors = find_neighbors_rect(idx_pair, n)
    idx = np.random.permutation(len(neighbors))[:np.minimum(option, len(neighbors))]
    for i in idx:
        neighbor = neighbors[i]
        if state[neighbor] == 0:
            state_new[neighbor] = 1

    return state_new


def update_Tm(state, idx_pair, n):  # migration
    state_new = state
    neighbors = find_neighbors_rect(idx_pair, n)
    # we use a random neighbor
    neighbor = neighbors[np.random.permutation(len(neighbors))[0]]
    if state[neighbor] == 0:
        state_new[idx_pair] = 0
        state_new[neighbor] = 1

    return state_new


def generate_from_probability(trans_1):
    ago = [np.sum(trans_1[:(i + 1)]) for i in range(len(trans_1))]
    ago.insert(0, 0)  # add 0 to generate roulette sequence
    ago = ago / np.sum(trans_1)  # normalization
    r = np.random.rand(1)  # uniform distribution
    for i in range(len(ago) - 1):
        if ago[i] <= r < ago[i + 1]:
            return i + 1  # return a chosen number (state)


def proliferation2inactive(state):
    i2x, i2y = np.where(state == 2)
    for x, y in zip(i2x, i2y):
        state[x, y] = 3
    return state


def migration2other(state, trans_1):
    i2x, i2y = np.where(state == 1)
    for x, y in zip(i2x, i2y):
        state[x, y] = generate_from_probability(trans_1)
    return state


def stimulation_v2(n: int, state_init, Tp1=2, Tp2=1, Tm=2, epoch=20, see_stage=False, trans_1=(2, 7, 1),
                          trans_2=(0, 1, 0)):
    # we call the overall cycle/period "epoch"
    # we call the three stages (proliferation 1,2; migration) "stage"
    # each stage involves Tp1/Tp2/Tm rounds
    # each round involves multiple "updates"

    states_epoch = np.zeros(shape=(epoch + 1, n, n), dtype='uint8')
    state_update = states_epoch[0] = state_init  # create a 3 dimensional array and initialize
    if see_stage:  # only if we want to observe the change after each stage, we add this, or it costs time
        states_stage = np.zeros(shape=(epoch * (Tp1 + Tp2 + Tm) + 1, n, n), dtype='uint8')
    else:
        states_stage = None

    for time in range(epoch):
        # state_init_round: keep the state at the start of this round, and use its 0/1/2
        # state_update: state after each update, continuously changing
        # copy means a no relationship with the original one

        for t in range(Tp1):
            state_init_round = state_update.copy()  # renew "init_round" with state_update each round
            idx_pairs = rand_pair(n)  # random order, no use!
            for idx_pair in idx_pairs:
                if state_init_round[idx_pair] == 2:  # proliferation state
                    state_update = update_Tp1(state_update, idx_pair, n, option=4)  # change "option" here manually
            if see_stage:
                states_stage[time * (Tp1 + Tp2 + Tm) + t + 1] = state_update.copy()

        for t in range(Tp2):
            state_init_round = state_update.copy()
            idx_pairs = rand_pair(n)
            for idx_pair in idx_pairs:
                if state_init_round[idx_pair] == 2:  # producing migration state
                    state_update = update_Tp2(state_update, idx_pair, n, option=4)
            if see_stage:
                states_stage[time * (Tp1 + Tp2 + Tm) + Tp1 + t + 1] = state_update.copy()

        proliferation2inactive(state_update)  # 2 -> 3

        for t in range(Tm):
            state_init_round = state_update.copy()
            idx_pairs = rand_pair(n)
            for idx_pair in idx_pairs:
                if state_init_round[idx_pair] == 1:  # migration state
                    state_update = update_Tm(state_update, idx_pair, n)
            if see_stage:
                states_stage[time * (Tp1 + Tp2 + Tm) + Tp1 + Tp2 + t + 1] = state_update.copy()

        state_update = migration2other(state_update, trans_1)  # 1 -> 1,2,3

        states_epoch[time + 1] = state_update
        if str(state_update).count("0") == 0:  # the squares are full; both the microorganisms
            print("the ground is full at the end of epoch {}".format(time + 1))
            return states_epoch[:time + 1], states_stage[:(time + 1) * (Tp1 + Tp2 + Tm)]

    return states_epoch, states_stage


# %% visualization
# see https://matplotlib.org/tutorials/colors/colormaps.html for more colormaps
def my_plot2d_animate(ca, title='evolve', interval=50, my_cmap='Greys'):
    cmap = plt.get_cmap(my_cmap)  # color
    fig = plt.figure()
    im = plt.imshow(ca[0], animated=True, cmap=cmap, vmin=0, vmax=3)  # set the value range here!
    plt.colorbar(ticks=[0, 1, 2, 3])  # set special values in the bar
    i = {'index': 0}

    def updatefig(*args):
        i['index'] += 1
        if i['index'] == len(ca):
            i['index'] = 0
        title = 'evolve at epoch ' + str(i['index'])  # manually set 'epoch' or 'stage'
        plt.title(title)
        im.set_array(ca[i['index']])
        return im,

    ani = animation.FuncAnimation(fig, updatefig, interval=interval, blit=True)
    # ani.save(filename="evolve.gif",writer='pillow')  # save as a gif file, change writer thus no error
    plt.show()


def my_plot2d(ca, timestep=None, title='', my_cmap='Greys'):
    cmap = plt.get_cmap(my_cmap)
    plt.title(title)
    if timestep is not None:
        data = ca[timestep]
    else:
        data = ca  # none represents a single state rather than a set of states
    plt.imshow(data, interpolation='none', cmap=cmap, vmin=0, vmax=3)
    plt.colorbar(ticks=[0, 1, 2, 3])
    plt.show()
