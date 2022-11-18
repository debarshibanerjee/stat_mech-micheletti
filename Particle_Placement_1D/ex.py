import numpy as np
from numpy.random import default_rng
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
from tqdm import tqdm


def get_distances(vec1, vec2, L):
    # dist_matrix = cdist(vec1[:, np.newaxis], vec2[:, np.newaxis], 'euclidean')
    dist_matrix = cdist(vec1[:, np.newaxis], vec2[:, np.newaxis], "minkowski", p=1)
    dist_matrix = np.abs(dist_matrix - L * (dist_matrix // (L / 2)))
    return dist_matrix


def generate_configs(num_trials, R_large, R_small, N_large, N_small, filling_fraction):

    L = ((2.0 * R_large * N_large) + (2.0 * R_small * N_small)) / filling_fraction

    rng = default_rng()
    macro_configs = rng.uniform(low=0.0, high=L, size=(num_trials, N_large))
    micro_configs = rng.uniform(low=0.0, high=L, size=(num_trials, N_small))

    results = []

    for i in tqdm(range(num_trials)):

        macro_vec = macro_configs[i, :]
        micro_vec = micro_configs[i, :]

        macro_macro_dist = get_distances(macro_vec, macro_vec, L)
        micro_micro_dist = get_distances(micro_vec, micro_vec, L)
        macro_micro_dist = get_distances(macro_vec, micro_vec, L)

        # print(macro_macro_dist)
        # print(macro_vec)
        # print(L)

        # Hack to ignore self-distances
        np.fill_diagonal(macro_macro_dist, np.max(macro_macro_dist))
        np.fill_diagonal(micro_micro_dist, np.max(micro_micro_dist))
        np.fill_diagonal(macro_micro_dist, np.max(macro_micro_dist))

        min_macro_macro_dist = np.min(macro_macro_dist)
        min_micro_micro_dist = np.min(micro_micro_dist)
        min_macro_micro_dist = np.min(macro_micro_dist)

        if (
            (min_macro_micro_dist < (2 * R_large))
            or (min_micro_micro_dist < (2 * R_small))
            or (min_macro_micro_dist < (R_large + R_small))
        ):
            continue
        else:
            # print(min_macro_macro_dist)
            results.append(min_macro_macro_dist)

    return results


def main():

    # parameters for macro particles
    N_large = 2
    R_large = 1.0

    # parameters for micro particles
    N_small = 20
    R_small = 0.01

    # phi
    filling_fraction = 0.3

    # Number of num_trials
    num_trials = 10000000

    # Number of runs
    num_runs = 10

    final_configs = []

    for i in range(num_runs):
        print(f"RUN {i}...")
        tmp_res = generate_configs(num_trials, R_large, R_small, N_large, N_small, filling_fraction)
        final_configs = final_configs + tmp_res

    print(f"Number of acceptable configurations = {len(final_configs)}")
    print(f"Total number of requested trials = {num_runs*num_trials}")

    np.savetxt("data.txt", final_configs)

    plt.hist(final_configs, bins="auto")
    plt.savefig("prob_distribution.png")
    plt.show()


## main() ##
if __name__ == "__main__":
    main()
