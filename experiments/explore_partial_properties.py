import numpy as np
import scanpy as sc
from src.paste.helper import kl_divergence, intersect, to_dense_array, extract_data_matrix
from src.paste.frantional_align import partial_pairwise_align
from helper import plot_slice
import matplotlib.pyplot as plt





def check_mapped_partial_region(sliceA_filename, sliceB_filename, m, alpha, armijo=False, dissimilarity='kl', norm=True):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])
    plot_slice(sliceA)
    plot_slice(sliceB)

    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))
    plot_slice(sliceA[common_spots])

    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=alpha, m=m, armijo=armijo, dissimilarity=dissimilarity, norm=norm, return_obj=True, verbose=True, matched_spots=matched_spots)
    print("Objective cost of the optimized alignment is %f" % log)
    print("Total mass transported is: " + str(np.sum(pi)))
    accuracy = 0
    for matched_spot in matched_spots:
        accuracy += pi[matched_spot[0]][matched_spot[1]]
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    print("Alignment accuracy is: " + str(accuracy))
    print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))


    going_out = np.dot(pi, np.ones((sliceB.shape[0], ))) > 0
    # going_out = np.sum(pi, axis=1) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in = np.dot(pi.T, np.ones((sliceA.shape[0], ))) > 0
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    plot_slice(going_out_part)
    plot_slice(coming_in_part)


    going_out_indices = np.argwhere(going_out > 0).flatten()
    coming_in_indices = np.argwhere(coming_in > 0).flatten()
    going_out_correct = 0
    coming_in_correct = 0
    for matched_spot in matched_spots:
        if matched_spot[0] in going_out_indices:
            going_out_correct += 1
        if matched_spot[1] in coming_in_indices:
            coming_in_correct += 1
    percentage_correctly_going_out = float(going_out_correct) / len(matched_spots)
    percentage_correctly_coming_in = float(coming_in_correct) / len(matched_spots)

    print("In the overlap region of the source slice, %f of them are actually mapped to some spot in the target slice" % percentage_correctly_going_out)
    print("In the overlap region of the target slice, %f of them are actually mapped to some spot in the source slice" % percentage_correctly_coming_in)

    plt.show()



#
# sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row0_col0.h5ad'
# sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row1_col0.h5ad'
#
# check_mapped_partial_region(sliceA_filename, sliceB_filename, m=0.5093484419263457, alpha=0.05, armijo=False, dissimilarity='kl', norm=True)
#
