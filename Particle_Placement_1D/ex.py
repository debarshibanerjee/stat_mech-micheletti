import concurrent.futures
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from numpy.random import default_rng
from scipy.spatial.distance import cdist
from tqdm import tqdm

# ================================================= Get Distances =====================================================#
def get_distances(vec1, vec2, L):
    """
    Function to return distance between 2 1-D arrays.
    Also takes care of PBC along 1-D, such that smaller distance is chosen.
    Returns the final distance matrix.
    """

    ## Get the distance matrix using either minkowski(p=1) or simple euclidean distance
    # dist_matrix = cdist(vec1[:, np.newaxis], vec2[:, np.newaxis], 'euclidean')
    dist_matrix = cdist(vec1[:, np.newaxis], vec2[:, np.newaxis], "minkowski", p=1)

    ## Take care of PBC: if dist > L/2, then choose the smaller distance
    ## // is floor division
    dist_matrix = np.abs(dist_matrix - L * (dist_matrix // (L / 2)))

    return dist_matrix


# =====================================================================================================================#


# ============================================ Generate Configurations ================================================#
def generate_configs(num_trials, R_large, R_small, N_large, N_small, filling_fraction):
    """
    Function to generate configurations, and keep those that meet the criteria.
    The rest are discarded.
    Returns the minimum distance between 2 large particles amongst the acceptable configurations.
    """

    ## Formula of L given as follows
    L = ((2.0 * R_large * N_large) + (2.0 * R_small * N_small)) / filling_fraction

    ## Recommended way of generating uniform random numbers acc to Numpy docs
    rng = default_rng()

    ## Generate appropriate configurations for both macro and micro particles
    ## Choose uniformly distributed random numbers between 0 and L
    ## Do so for all possible trials
    macro_configs = rng.uniform(low=0.0, high=L, size=(num_trials, N_large))
    micro_configs = rng.uniform(low=0.0, high=L, size=(num_trials, N_small))

    ## List to append results
    results = []

    ## Start the loop over all trials
    for i in tqdm(range(num_trials)):

        ## Choose the configurations for a given trial
        macro_vec = macro_configs[i, :]
        micro_vec = micro_configs[i, :]

        ## Get the distance matrices for macro-macro, micro-micro, and macro-micro combinations
        macro_macro_dist = get_distances(macro_vec, macro_vec, L)
        micro_micro_dist = get_distances(micro_vec, micro_vec, L)
        macro_micro_dist = get_distances(macro_vec, micro_vec, L)

        ## Some debugging comments
        # print(macro_macro_dist)
        # print(macro_vec)
        # print(L)

        ## Hack to ignore self-distances by setting diagonal elements to largest value
        ## This is because diagonal elements are otherwise set to 0 for obvious reasons
        ## Namely, that they are distances between 1 particle and itself
        np.fill_diagonal(macro_macro_dist, np.max(macro_macro_dist))
        np.fill_diagonal(micro_micro_dist, np.max(micro_micro_dist))
        np.fill_diagonal(macro_micro_dist, np.max(macro_micro_dist))

        ## Some debugging comments
        # print(macro_macro_dist)
        # print(np.unique(macro_macro_dist))

        ## Choose minimum distances for all cases of particles
        min_macro_macro_dist = np.min(macro_macro_dist)
        min_micro_micro_dist = np.min(micro_micro_dist)
        min_macro_micro_dist = np.min(macro_micro_dist)

        ## If the spheres oriented along 1D overlap each other, ignore the configurations
        ## In that case, move on to the next loop iteration, i.e., the next trial
        ## The distance must be less than sum of radii of 2 compared spheres for configuration to be unacceptable
        ## If all the checks pass, add the minimum distance between 2 large particles to a List
        ## Remember that ultimately we are interested in the probability distribution of this distance for macro particles
        if (
            (min_macro_macro_dist < (2 * R_large))
            or (min_micro_micro_dist < (2 * R_small))
            or (min_macro_micro_dist < (R_large + R_small))
        ):
            continue
        else:
            # print(min_macro_macro_dist)
            results.append(np.unique(macro_macro_dist))

    return results


# =====================================================================================================================#


# ================================================= Define main() =====================================================#
def main():

    n_args = len(sys.argv)

    if n_args != 2:
        nprocs = 1
        print("You can also run the code as: 'python3 ex.py 4', to use 4 CPUs")
    else:
        nprocs = int(sys.argv[1])

    n_cpus = os.cpu_count()
    if nprocs > n_cpus:
        nprocs = n_cpus

    print(f"Number of processes created simultaneously = {nprocs}")

    ## parameters for macro particles
    N_large = 2
    R_large = 1.0

    ## parameters for micro particles
    N_small = 20
    R_small = 0.01

    ## phi
    filling_fraction = 0.5

    ## Number of num_trials
    num_trials = 1000000

    ## Number of runs
    num_runs = 300

    tmp_configs = []
    final_configs = []

    ## Parallel execution
    futures = []
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=nprocs)
    for i in range(num_runs):
        print(f"RUN {i+1}...")
        futures.append(
            executor.submit(generate_configs, num_trials, R_large, R_small, N_large, N_small, filling_fraction)
        )

    ## Recover results from futures
    for future in concurrent.futures.as_completed(futures):
        tmp_configs.append(future.result())

    for i in range(len(tmp_configs)):
        l = [x for x in tmp_configs[i]]
        final_configs = final_configs + l

    print(f"Number of acceptable configurations = {len(final_configs)}")
    print(f"Total number of requested trials = {num_runs*num_trials}")

    np.savetxt("data.txt", final_configs)

    plt.hist(final_configs, bins="auto", density=True)
    plt.savefig("prob_distribution.png")
    # plt.show()


# =====================================================================================================================#


# ================================================== Call main() ======================================================#
if __name__ == "__main__":
    main()
# =====================================================================================================================#
