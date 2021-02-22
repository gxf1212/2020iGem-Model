# growth, to spores, nutrient consumption, rps production (with N), carrying capacity
def v0_10_22():
    probs_migration = [0.02, 0.002, 0.002, 0.02]  # diffusion coefficient
    weight = [[1, 0.05, 1], [1, 0.05, 1]]
    params_BS = (0.40, 0.057, 0.1, 0.1, 1000)  # ?
    params_No = (0.01, 0.001, -0.1, 0.5, 1000)
    init_classic = (100, 200, 200, 1000)
    return probs_migration, weight, params_BS, params_No, init_classic


def v1_high_no():
    probs_migration = [0.02, 0.002, 0.002, 0.02]  # diffusion coefficient
    weight = [[1, 0.05, 1], [1, 0.05, 1]]
    params_BS = (0.35, 0.05, 0.1, 0.1, 1000)
    params_No = (0.3, 0.03, -0.2, 0.5, 1000)
    # params_No = (0.3, 0.03, -0.2, 0.5, 1000)
    init_classic = (100, 200, 200, 1000)
    return probs_migration, weight, params_BS, params_No, init_classic


def v2_low_nut():
    probs_migration = [0.02, 0.008, 0.002, 0.02]
    weight = [[1, 0.05, 1], [1, 0.05, 1]]
    params_BS = (0.35, 0.05, 1, 1, 1000)
    params_No = (0.3, 0.03, -2, 5, 1000)
    init_classic = (100, 200, 200, 10000)
    return probs_migration, weight, params_BS, params_No, init_classic


def v3_to_nut():
    probs_migration = [0.02, 0.002, 0.002, 0.02]  # diffusion coefficient
    weight = [[1, 0.05, -0.2], [1, 0.05, -0.2]]
    params_BS = (0.35, 0.05, 0.1, 0.1, 1000)
    params_No = (0.3, 0.03, -0.2, 0.5, 1000)
    # params_No = (0.3, 0.03, -0.2, 0.5, 1000)
    init_classic = (100, 200, 200, 1000)
    return probs_migration, weight, params_BS, params_No, init_classic


def v3_low_nut():
    probs_migration = [0.02, 0.008, 0.002, 0.02]
    weight = [[1, 0.05, 1], [1, 0.05, 1]]
    params_BS = (0.35, 0.05, 1, 0.1, 1000)
    params_No = (0.3, 0.03, -2, 0.5, 1000)
    init_classic = (100, 200, 200, 1000)
    return probs_migration, weight, params_BS, params_No, init_classic


def v4_high_rps():
    probs_migration = [0.02, 0.002, 0.002, 0.02]  # diffusion coefficient
    weight = [[1, 0.05, 1], [1, 0.05, 1]]
    params_BS = (0.35, 0.05, 0.1, 0.5, 1000)
    params_No = (0.3, 0.03, -0.2, 2.5, 1000)
    # params_No = (0.3, 0.03, -0.2, 0.5, 1000)
    init_classic = (100, 200, 200, 1000)
    return probs_migration, weight, params_BS, params_No, init_classic


def vn_real():
    # real params
    # these params should be proportional to the original
    # growth/decay rate: *20
    # params_BS = (np.log(2)/30*15, 0.0033*15, 1, 0.1, 1000)
    # params_No = (6.5917e-04*15, 5.83e-05*15, -2, 0.1, 1000)
    # params_BS = (np.log(2)/30*15, 0.0033*15, 1, 0.1, 1000)
    # params_No = (6.5917e-04*225, 5.83e-05*225, -2, 0.1, 1000)

    return


def vn_try():
    params_BS = (0.35, 0.05, 1, 0.1, 5000)
    params_No = (0.3, 0.03, -2, 0.1, 5000)  # keep r/m ratio
    return

