import numpy as np
import scanpy as sc
from src.paste.helper import kl_divergence, intersect, to_dense_array, extract_data_matrix
from src.paste.fractional_align import partial_pairwise_align, partial_pairwise_align_paste_init
from experiments.helper import plot_slice, plot_slice_umi, compute_alignment_ari, plot_slices_overlap
import matplotlib.pyplot as plt
import src.paste.PASTE as paste
from src.paste.visualization import partial_stack_slices_pairwise, stack_slices_pairwise




def check_mapped_partial_region(sliceA_filename, sliceB_filename, m, alpha, armijo=False, dissimilarity='kl', norm=True):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])
    plot_slice(sliceA)
    plot_slice_umi(sliceA)
    plot_slice(sliceB)
    plot_slice_umi(sliceB)

    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))
    # plot_slice(sliceA[common_spots])
    # plot_slice_umi(sliceA[common_spots])


    # pi, log = partial_pairwise_align(sliceA, sliceB, alpha=alpha, m=m, armijo=armijo, dissimilarity=dissimilarity, norm=norm, return_obj=True, verbose=True, matched_spots=matched_spots)
    pi, log = partial_pairwise_align_paste_init(sliceA, sliceB, alpha=alpha, m=m, armijo=armijo, dissimilarity=dissimilarity, norm=norm, return_obj=True, verbose=True, matched_spots=matched_spots)
    print(pi.shape)
    print("Objective cost of the optimized alignment is %f" % log)
    print("Total mass transported is: " + str(np.sum(pi)))
    accuracy = 0
    for matched_spot in matched_spots:
        accuracy += pi[matched_spot[0]][matched_spot[1]]
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    print("Alignment accuracy is: " + str(accuracy))
    print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))
    ari = compute_alignment_ari(sliceA, sliceB, pi)
    print("ARI is: " + str(ari))


    going_out = np.dot(pi, np.ones((sliceB.shape[0], ))) > 0
    # going_out = np.sum(pi, axis=1) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in = np.dot(pi.T, np.ones((sliceA.shape[0], ))) > 0
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    plot_slice(going_out_part)
    # plot_slice_umi(going_out_part)
    plot_slice(coming_in_part)
    # plot_slice_umi(coming_in_part)


    # going_out_indices = np.argwhere(going_out > 0).flatten()
    # coming_in_indices = np.argwhere(coming_in > 0).flatten()
    # going_out_correct = 0
    # coming_in_correct = 0
    # for matched_spot in matched_spots:
    #     if matched_spot[0] in going_out_indices:
    #         going_out_correct += 1
    #     if matched_spot[1] in coming_in_indices:
    #         coming_in_correct += 1
    # percentage_correctly_going_out = float(going_out_correct) / len(matched_spots)
    # percentage_correctly_coming_in = float(coming_in_correct) / len(matched_spots)
    # print("In the overlap region of the source slice, %f of them are actually mapped to some spot in the target slice" % percentage_correctly_going_out)
    # print("In the overlap region of the target slice, %f of them are actually mapped to some spot in the source slice" % percentage_correctly_coming_in)


    # overlap_source_indices = [matched_spot[0] for matched_spot in matched_spots]
    # overlap_dest_indices = [matched_spot[1] for matched_spot in matched_spots]
    # cnt_source_off = 0
    # cnt_target_off = 0
    # for source_spot_idx in overlap_source_indices:
    #     mapped_to_indices = np.argwhere(pi[source_spot_idx] > 0).flatten()
    #     for idx in mapped_to_indices:
    #         if idx not in overlap_dest_indices:
    #             cnt_source_off += 1
    #             break
    # percentage_source_off = float(cnt_source_off) / len(overlap_source_indices)
    # for dest_spot_idx in overlap_dest_indices:
    #     mapped_to_indices = np.argwhere(pi[:, dest_spot_idx] > 0).flatten()
    #     for idx in mapped_to_indices:
    #         if idx not in overlap_source_indices:
    #             cnt_target_off += 1
    #             break
    # percentage_target_off = float(cnt_target_off) / len(overlap_dest_indices)
    # print("In the overlap region of the source slice, %f of them have mass transported to spots in the un-overlapped region in the target slice" % percentage_source_off)
    # print("In the overlap region of the target slice, %f of them have mass transported to spots in the un-overlapped region in the source slice" % percentage_target_off)


    # umi_of_correct_spots = []
    # umi_of_wrong_spots = []
    # for matched_spot in matched_spots:
    #     if pi[matched_spot[0]][matched_spot[1]] > 0:
    #         umi_of_correct_spots.append(sliceA.obs["sum_umi"][matched_spot[0]])
    #     else:
    #         umi_of_wrong_spots.append(sliceA.obs["sum_umi"][matched_spot[0]])
    # avg_umi_of_correctly_mapped_spots = np.average(umi_of_correct_spots)
    # avg_umi_of_wrongly_unmapped_spots = np.average(umi_of_wrong_spots)
    # print("In the overlap region of the source slice, the average UMI of the correctly mapped spots is %f" % avg_umi_of_correctly_mapped_spots)
    # print("In the overlap region of the source slice, the average UMI of unmapped spots is %f" % avg_umi_of_wrongly_unmapped_spots)

    source_correctly_mapped = []
    target_correctly_mapped = []
    for matched_spot in matched_spots:
        if pi[matched_spot[0]][matched_spot[1]] > 0:
            source_correctly_mapped.append(matched_spot[0])
            target_correctly_mapped.append(matched_spot[1])
    source_correctly_mapped_part = sliceA[sliceA.obs.index[np.sort(source_correctly_mapped)]]
    target_correctly_mapped_part = sliceB[sliceB.obs.index[np.sort(target_correctly_mapped)]]
    plot_slice(source_correctly_mapped_part)
    plot_slice(target_correctly_mapped_part)


    plt.show()


#
# sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row0_col0.h5ad'
# sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row1_col0.h5ad'
#
# check_mapped_partial_region(sliceA_filename, sliceB_filename, m=0.5093484419263457, alpha=0.05, armijo=False, dissimilarity='kl', norm=True)
#


def original_paste_inspect_pi(sliceA_filename, sliceB_filename, alpha, dissimilarity):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])

    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))

    pi, log = paste.pairwise_align(sliceA, sliceB, alpha=alpha, dissimilarity=dissimilarity, norm=True, return_obj=True,
                                   verbose=True)
    print(pi.shape)

    # Calculate accuracy and maximum possible accuracy
    accuracy = 0
    for matched_spot in matched_spots:
        accuracy += pi[matched_spot[0]][matched_spot[1]]
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    print("Alignment accuracy is: " + str(accuracy))
    print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))

    # Get distribution of number of spots in destination that one spot in source is mapped to, for spots in overlapped
    # and non-overlapped region
    overlap_source_indices = [matched_spot[0] for matched_spot in matched_spots]
    overlap_distribution = []
    nonoverlap_distribution = []
    for source_spot_idx in overlap_source_indices:
        going_out = pi[source_spot_idx]
        overlap_distribution.append((going_out > 0).sum())
    for i in range(pi.shape[0]):
        if i not in overlap_source_indices:
            going_out = pi[i]
            nonoverlap_distribution.append((going_out > 0).sum())

    return overlap_distribution, nonoverlap_distribution


# sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row0_col0.h5ad'
# sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row1_col0.h5ad'
#
# overlap_distribution, nonoverlap_distribution = original_paste_inspect_pi(sliceA_filename, sliceB_filename, alpha=0.05, dissimilarity='kl')
# print(overlap_distribution)
# print(nonoverlap_distribution)
# print(len(overlap_distribution))
# print(len(nonoverlap_distribution))





def check_mapped_partial_region_realdata(sliceA_filename, sliceB_filename, m, alpha, armijo=False, dissimilarity='glmpca', norm=True):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])
    plot_slice(sliceA)
    plot_slice_umi(sliceA)
    plot_slice(sliceB)
    plot_slice_umi(sliceB)

    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=alpha, m=m, armijo=armijo, dissimilarity=dissimilarity, norm=norm, return_obj=True, verbose=True)
    print(pi.shape)
    print("Total mass transported is: " + str(np.sum(pi)))
    print("ARI is: " + str(compute_alignment_ari(sliceA, sliceB, pi)))

    going_out = np.dot(pi, np.ones((sliceB.shape[0], ))) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in = np.dot(pi.T, np.ones((sliceA.shape[0], ))) > 0
    # print(np.sum(np.dot(pi.T, np.ones((sliceA.shape[0], )))))
    # print(np.max(np.dot(pi.T, np.ones((sliceA.shape[0],)))))
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    plot_slice(going_out_part)
    plot_slice(coming_in_part)

    """
    stacking
    """
    new_slices = partial_stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()



def original_realdata(sliceA_filename, sliceB_filename, alpha, dissimilarity='kl', norm=True):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    plot_slice(sliceA)
    plot_slice(sliceB)
    pi, log = paste.pairwise_align(sliceA, sliceB, alpha=alpha, dissimilarity=dissimilarity, norm=norm, return_obj=True, verbose=True)
    print("ARI is: " + str(compute_alignment_ari(sliceA, sliceB, pi)))

    """
    stacking
    """
    new_slices = stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()



