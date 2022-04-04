import numpy as np
import scanpy as sc
from src.paste.fractional_align import partial_pairwise_align
from src.paste.helper import intersect
from experiments.helper import compute_alignment_ari



def partial_paste_pairwise_align(sliceA_filename, sliceB_filename, m, alpha, armijo=False, dissimilarity='glmpca', norm=True):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])

    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))

    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=alpha, m=m, armijo=armijo, dissimilarity=dissimilarity, norm=norm, return_obj=True, verbose=True)
    # print("Objective cost is: " + str(log))
    print("Total mass transported is: " + str(np.sum(pi)))
    ari = compute_alignment_ari(sliceA, sliceB, pi)

    accuracy = 0
    for matched_spot in matched_spots:
        accuracy += pi[matched_spot[0]][matched_spot[1]]
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    return maximum_possible_accuracy, accuracy, ari, log



m_to_test = [0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/pca/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta0.1_row0_col0.h5ad'
sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/pca/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta0.1_row1_col0.h5ad'
for m in m_to_test:
    print("=======================================")
    maximum_possible_accuracy, accuracy, ari, log = partial_paste_pairwise_align(sliceA_filename, sliceB_filename, m, alpha=0.1)
    print("m = " + str(m))
    print("Maximum possible accuracy = " + str(maximum_possible_accuracy))
    print("Accuracy = " + str(accuracy))
    print("ARI = " + str(ari))
    print("Objective cost = " + str(log))

