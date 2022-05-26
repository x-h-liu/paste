import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from src.paste.fractional_align import partial_pairwise_align_given_cost_matrix
from src.paste.helper import intersect, to_dense_array, extract_data_matrix, glmpca_distance
from experiments.helper import compute_alignment_ari


# m_to_run = [0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
m_to_run = [0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]


def get_cost_curve(m_to_run):
    #sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/pca/151674_overlap0.9_dropFalse_rotateFalse_resampleTrue_delta1.0_row0_col0.h5ad'
    #sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/pca/151674_overlap0.9_dropFalse_rotateFalse_resampleTrue_delta1.0_row1_col0.h5ad'
    sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151673/151673_overlap0.5_dropFalse_rotateFalse_resampleFalse_delta0.0_row0_col0.h5ad'
    sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674/151674_overlap0.5_dropFalse_rotateFalse_resampleFalse_delta0.0_row1_col0.h5ad'
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    common_genes = intersect(sliceA.var.index, sliceB.var.index)
    sliceA = sliceA[:, common_genes]
    sliceB = sliceB[:, common_genes]
    # Get transport cost matrix
    A_X, B_X = to_dense_array(extract_data_matrix(sliceA, None)), to_dense_array(extract_data_matrix(sliceB, None))
    M = glmpca_distance(A_X, B_X, latent_dim=50, filter=True)

    # Get variables for performance evaluation
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])
    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)

    costs = []
    for m in m_to_run:
        print("============================")
        pi, log = partial_pairwise_align_given_cost_matrix(sliceA, sliceB, M=M, alpha=0.5, m=m, armijo=False,
                                                           norm=True, return_obj=True, verbose=True)
        costs.append(log)
        ari = compute_alignment_ari(sliceA, sliceB, pi)
        accuracy = 0
        for matched_spot in matched_spots:
            accuracy += pi[matched_spot[0]][matched_spot[1]]
        print("m = " + str(m))
        print("Total mass transported is: " + str(np.sum(pi)))
        print("Maximum possible accuracy = " + str(maximum_possible_accuracy))
        print("Accuracy = " + str(accuracy))
        print("ARI = " + str(ari))
        print("Objective cost = " + str(log))
    print("============================")
    print(costs)
    return costs


def plot_cost_curve(m_list, cost_list):
    fig, ax = plt.subplots()
    ax.plot(m_list, cost_list)
    ax.set_xlim(1, 0)
    ax.set_xticks([0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
    ax.set_xlabel("m")
    ax.set_ylabel("Objective cost")
    # ax.set_title("Delta=1.0, Truth=0.9")
    ax.set_title("real0001, alpha=0, Truth=0.5")
    plt.show()


get_cost_curve(m_to_run)
# cost_curve = [17.594347785295202, 16.47327549096224, 15.235313195932259, 14.090370372643813, 13.015117539229147, 11.991864219853296, 11.010487609492362, 10.066521592246898, 9.154409734954998, 8.267415934498535, 7.40681850288412, 6.568295547351732, 5.749842361970782, 4.949811582112172, 4.170637869943312, 3.411097723554066, 2.6727839708318686, 1.9565989940666406, 1.2659754390320561, 0.6053882741180994]
#
# plot_cost_curve(m_to_run, cost_curve)

