import numpy as np
import scanpy as sc
import anndata as ad
import src.paste.PASTE as paste
from src.paste.helper import kl_divergence, intersect, to_dense_array, extract_data_matrix
from src.paste.frantional_align import partial_pairwise_align, partial_pairwise_align_POT


#sliceA = sc.read_h5ad('/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row0_col0.h5ad')
# print(sliceA)
# print(sliceA.X)
# print(sliceA.var.index)
# print(sliceA.obsm['spatial'])
# print(sliceA.obsm['spatial'].shape)
# spotnames = sliceA.obs.index
# print(spotnames)
# datamatrix = to_dense_array(sliceA.X)
#
# x = datamatrix[spotnames.get_loc('TTGTTGTGTGTCAAGA-1.9')]
# print(x)
# print("===========")
# y = to_dense_array(sliceA['TTGTTGTGTGTCAAGA-1.9'].X)[0]
# print(y)
# print(np.array_equal(x, y))


# print(spotnames.get_loc('AAACAAGTATCTCCCA-1.9'))
# print(spotnames.get_loc('AAACAATCTACTAGCA-1.4'))
# print(spotnames.get_loc('TTGTTTGTGTAAATTC-1.9'))

#sliceB = sc.read_h5ad('/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row0_col1.h5ad')
# print(sliceB)
# common_genes = intersect(sliceA.var.index, sliceB.var.index)
# print(type(common_genes))
# print(common_genes)
# print(sliceA[:, common_genes])
# print(sliceA[:, common_genes].shape)
# print(sliceB.X.shape)

# maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])
# spotnamesA = sliceA.obs.index
# spotnamesB = sliceB.obs.index
# common_spots = intersect(spotnamesA, spotnamesB)
# matched_spots = []
# for spot in common_spots:
#     matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))
#
# pi, log = paste.pairwise_align(sliceA, sliceB, return_obj=True)
# print(pi.shape)
# print(log)
#
# accuracy = 0
# for matched_spot in matched_spots:
#     accuracy += pi[matched_spot[0]][matched_spot[1]]
#
# maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
#
# print(accuracy)
# print(maximum_possible_accuracy)


def original_paste_pairwise_align(sliceA_filename, sliceB_filename, alpha, dissimilarity):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])

    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))

    pi, log = paste.pairwise_align(sliceA, sliceB, alpha=alpha, dissimilarity=dissimilarity, norm=True, return_obj=True, verbose=True)
    # print("Alignment value is: " + str(log))
    # for loss in log['loss']:
    #     print(loss)

    accuracy = 0
    for matched_spot in matched_spots:
        accuracy += pi[matched_spot[0]][matched_spot[1]]
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    return accuracy, maximum_possible_accuracy


# accuracy, maximum_possible_accuracy = original_paste_pairwise_align('/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row0_col1.h5ad',
#                               '/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row1_col1.h5ad')
#
# print(accuracy, maximum_possible_accuracy)


# to_align = [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (1, 1))]
# #to_align = [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (1, 1)), ((0, 0), (1, 1)), ((0, 1), (1, 0))]
# #to_align = [((0, 0), (1, 0))]
# for pair in to_align:
#     sliceA_row = pair[0][0]
#     sliceA_col = pair[0][1]
#     sliceB_row = pair[1][0]
#     sliceB_col = pair[1][1]
#     # sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Data/PASTE/delta1/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_delta1_row' \
#     #                   + str(sliceA_row) + '_col' + str(sliceA_col) + '.h5ad'
#     # sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Data/PASTE/delta1/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_delta1_row' \
#     #                   + str(sliceB_row) + '_col' + str(sliceB_col) + '.h5ad'
#     sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674_overlap0.3_dropFalse_rotateFalse_resampleTrue_delta0.5_row' \
#                       + str(sliceA_row) + '_col' + str(sliceA_col) + '.h5ad'
#     sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674_overlap0.3_dropFalse_rotateFalse_resampleTrue_delta0.5_row' \
#                       + str(sliceB_row) + '_col' + str(sliceB_col) + '.h5ad'
#     print("=======================")
#     print(pair)
#     accuracy, maximum_possible_accuracy = original_paste_pairwise_align(sliceA_filename, sliceB_filename, alpha=0.05, dissimilarity='kl')
#     print("Alignment accuracy is: " + str(accuracy))
#     print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))


"""
Code for bulk running Original
"""
# overlap_to_run = [0.9, 0.7, 0.5, 0.3]
# # delta_to_run = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0]
# delta_to_run = [1.0]
# # alpha_to_run = [0.0, 0.05, 0.1, 0.2]
# alpha_to_run = [0.1]
# to_align = [((0, 0), (1, 0))]
#
# overlap_alpha_to_delta_accuracy = {}
# for overlap in overlap_to_run:
#     overlap_alpha_to_delta_accuracy[overlap] = {}
#     for alpha in alpha_to_run:
#         overlap_alpha_to_delta_accuracy[overlap][alpha] = []
#         for delta in delta_to_run:
#             for pair in to_align:
#                 sliceA_row = pair[0][0]
#                 sliceA_col = pair[0][1]
#                 sliceB_row = pair[1][0]
#                 sliceB_col = pair[1][1]
#                 sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151510/151510_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row' \
#                                   + str(sliceA_row) + '_col' + str(sliceA_col) + '.h5ad'
#                 sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151510/151510_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row' \
#                                   + str(sliceB_row) + '_col' + str(sliceB_col) + '.h5ad'
#                 print("=======================")
#                 # print(pair)
#                 print("Overlap={0}, Alpha={1}, Delta={2}".format(overlap, alpha, delta))
#                 accuracy, maximum_possible_accuracy = original_paste_pairwise_align(sliceA_filename, sliceB_filename,
#                                                                                     alpha=alpha, dissimilarity='kl')
#                 print("Alignment accuracy is: " + str(accuracy))
#                 print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))
#                 overlap_alpha_to_delta_accuracy[overlap][alpha].append(accuracy)
#
#
# print("\n****************************************************************************")
# for overlap in overlap_alpha_to_delta_accuracy:
#     for alpha in overlap_alpha_to_delta_accuracy[overlap]:
#         print(overlap, alpha)
#         print(overlap_alpha_to_delta_accuracy[overlap][alpha])
"""
Code for bulk running Original ends
"""


def partial_POT_pairwise_align(sliceA_filename, sliceB_filename, m, dissimilarity='kl', norm=False):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])

    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))

    pi, log = partial_pairwise_align_POT(sliceA, sliceB, m=m, dissimilarity=dissimilarity, norm=norm, return_obj=True)
    print("Alignment value is: " + str(log))
    print("Total mass transported is: " + str(np.sum(pi)))

    accuracy = 0
    for matched_spot in matched_spots:
        accuracy += pi[matched_spot[0]][matched_spot[1]]
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    return accuracy, maximum_possible_accuracy



def partial_paste_pairwise_align(sliceA_filename, sliceB_filename, m, alpha, armijo=False, dissimilarity='kl', norm=False):
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
    # print("Alignment value is: " + str(log))
    # for loss in log['loss']:
    #     print(loss)
    print("Total mass transported is: " + str(np.sum(pi)))

    accuracy = 0
    for matched_spot in matched_spots:
        accuracy += pi[matched_spot[0]][matched_spot[1]]
    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    return accuracy, maximum_possible_accuracy

# accuracy, maximum_possible_accuracy = partial_paste_pairwise_align('/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_row0_col0.h5ad',
#                               '/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_row1_col0.h5ad',
#                                                                    m=0.8, alpha=0.1, armijo=False)
#
# print(accuracy, maximum_possible_accuracy)


# # to_align = [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (1, 1)), ((0, 0), (1, 1)), ((0, 1), (1, 0))]
# # to_align = [((0, 0), (0, 1)), ((0, 0), (1, 0)), ((0, 1), (1, 1)), ((1, 0), (1, 1))]
# to_align = [((0, 0), (1, 0))]
# for pair in to_align:
#     sliceA_row = pair[0][0]
#     sliceA_col = pair[0][1]
#     sliceB_row = pair[1][0]
#     sliceB_col = pair[1][1]
#     # sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Data/PASTE/delta1/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_delta1_row' \
#     #                   + str(sliceA_row) + '_col' + str(sliceA_col) + '.h5ad'
#     # sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Data/PASTE/delta1/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_delta1_row' \
#     #                   + str(sliceB_row) + '_col' + str(sliceB_col) + '.h5ad'
#     sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674_overlap0.7_dropFalse_rotateFalse_resampleTrue_delta0.5_row' \
#                       + str(sliceA_row) + '_col' + str(sliceA_col) + '.h5ad'
#     sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674_overlap0.7_dropFalse_rotateFalse_resampleTrue_delta0.5_row' \
#                       + str(sliceB_row) + '_col' + str(sliceB_col) + '.h5ad'
#     print("=======================")
#     print(pair)
#     accuracy, maximum_possible_accuracy = partial_paste_pairwise_align(sliceA_filename, sliceB_filename, m=0.75, alpha=0.05,
#                                                                        armijo=False, dissimilarity='kl', norm=True)
#     print("Alignment accuracy is: " + str(accuracy))
#     print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))


"""
Code for bulk running Partial
"""
overlap_to_run = [0.9, 0.7, 0.5, 0.3]
# delta_to_run = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0]
delta_to_run = [0.0]
# alpha_to_run = [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
alpha_to_run = [0.05]
to_align = [((0, 0), (1, 0))]


overlap_alpha_to_delta_accuracy = {}
for overlap in overlap_to_run:
    overlap_alpha_to_delta_accuracy[overlap] = {}
    for alpha in alpha_to_run:
        overlap_alpha_to_delta_accuracy[overlap][alpha] = []
        for delta in delta_to_run:
            for pair in to_align:
                sliceA_row = pair[0][0]
                sliceA_col = pair[0][1]
                sliceB_row = pair[1][0]
                sliceB_col = pair[1][1]
                sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row' \
                                  + str(sliceA_row) + '_col' + str(sliceA_col) + '.h5ad'
                sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row' \
                                  + str(sliceB_row) + '_col' + str(sliceB_col) + '.h5ad'
                print("=======================")
                # print(pair)
                print("Overlap={0}, Alpha={1}, Delta={2}".format(overlap, alpha, delta))
                accuracy, maximum_possible_accuracy = partial_paste_pairwise_align(sliceA_filename, sliceB_filename, m=overlap, alpha=alpha, armijo=False, dissimilarity='kl', norm=True)
                print("Alignment accuracy is: " + str(accuracy))
                print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))
                overlap_alpha_to_delta_accuracy[overlap][alpha].append(accuracy)


print("\n****************************************************************************")
for overlap in overlap_alpha_to_delta_accuracy:
    for alpha in overlap_alpha_to_delta_accuracy[overlap]:
        print(overlap, alpha)
        print(overlap_alpha_to_delta_accuracy[overlap][alpha])
"""
Code for bulk running Partial ends
"""



"""
Code for bulk running Partial POT
"""
# overlap_to_run = [0.9, 0.7, 0.5, 0.3]
# delta_to_run = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0]
# alpha_to_run = [0.0]
# to_align = [((0, 0), (1, 0))]
#
# print("PARTIAL POT")
# overlap_alpha_to_delta_accuracy = {}
# for overlap in overlap_to_run:
#     overlap_alpha_to_delta_accuracy[overlap] = {}
#     for alpha in alpha_to_run:
#         overlap_alpha_to_delta_accuracy[overlap][alpha] = []
#         for delta in delta_to_run:
#             for pair in to_align:
#                 sliceA_row = pair[0][0]
#                 sliceA_col = pair[0][1]
#                 sliceB_row = pair[1][0]
#                 sliceB_col = pair[1][1]
#                 sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row' \
#                                   + str(sliceA_row) + '_col' + str(sliceA_col) + '.h5ad'
#                 sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row' \
#                                   + str(sliceB_row) + '_col' + str(sliceB_col) + '.h5ad'
#                 print("=======================")
#                 print("Overlap={0}, Alpha={1}, Delta={2}".format(overlap, alpha, delta))
#                 accuracy, maximum_possible_accuracy = partial_POT_pairwise_align(sliceA_filename, sliceB_filename, m=overlap, dissimilarity='kl', norm=True)
#                 print("Alignment accuracy is: " + str(accuracy))
#                 print("Maximum possible accuracy is: " + str(maximum_possible_accuracy))
#                 overlap_alpha_to_delta_accuracy[overlap][alpha].append(accuracy)
#
#
# print("\n****************************************************************************")
# for overlap in overlap_alpha_to_delta_accuracy:
#     for alpha in overlap_alpha_to_delta_accuracy[overlap]:
#         print(overlap, alpha)
#         print(overlap_alpha_to_delta_accuracy[overlap][alpha])
"""
Code for bulk running Partial POT ends
"""