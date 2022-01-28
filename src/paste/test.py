import numpy as np
import scanpy as sc
import anndata as ad
import src.paste.PASTE as paste
from src.paste.helper import kl_divergence, intersect, to_dense_array, extract_data_matrix
from scipy.spatial import distance_matrix


# def gwgrad_partial(C1, C2, T):
#     """Compute the GW gradient. Note: we can not use the trick in :ref:`[12] <references-gwgrad-partial>`
#     as the marginals may not sum to 1.
#
#     Parameters
#     ----------
#     C1: array of shape (n_p,n_p)
#         intra-source (P) cost matrix
#
#     C2: array of shape (n_u,n_u)
#         intra-target (U) cost matrix
#
#     T : array of shape(n_p+nb_dummies, n_u) (default: None)
#         Transport matrix
#
#     Returns
#     -------
#     numpy.array of shape (n_p+nb_dummies, n_u)
#         gradient
#
#
#     .. _references-gwgrad-partial:
#     References
#     ----------
#     .. [12] Peyré, Gabriel, Marco Cuturi, and Justin Solomon,
#         "Gromov-Wasserstein averaging of kernel and distance matrices."
#         International Conference on Machine Learning (ICML). 2016.
#     """
#     cC1 = np.dot(C1 ** 2 / 2, np.dot(T, np.ones(C2.shape[0]).reshape(-1, 1)))
#     print(cC1)
#     cC2 = np.dot(np.dot(np.ones(C1.shape[0]).reshape(1, -1), T), C2 ** 2 / 2)
#     print(cC2)
#     constC = cC1 + cC2
#     print(constC)
#
#     A = np.dot(cC1, np.ones((1, C2.shape[0])))
#     print(A)
#     B = np.dot(np.ones((C1.shape[0], 1)), cC2)
#     print(B)
#     print(A + B)
#
#
#
#
#
#     A = -np.dot(C1, T).dot(C2.T)
#     tens = constC + A
#     return tens * 2
#
#
# p = 2
# q = 3
#
# C1 = np.ones((p, p))
# C2 = np.ones((q, q))
# T = np.ones((p, q))
# gwgrad_partial(C1, C2, T)
#
#
#
#
# A = np.array([[2,2],[2,2]])
# print(A)
# print(A**3)


# sliceA = sc.read_10x_h5("/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row0_col0.h5ad")

#sliceA = ad.read('/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row0_col0.h5ad')
# print(sliceA)
# print(sliceA.var.index)
# print(sliceA.obs)
# print(sliceA.obs.index)
# print(sliceA.obsm['spatial'])
# print(sliceA.obsm['spatial'].shape)
# print(sliceA.X)


# sliceB = ad.read('/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row0_col1.h5ad')
# print(sliceB)
#
#
# # sliceC = ad.read('/Users/xinhaoliu/Desktop/Research/Data/PASTE/Share/151674_overlap1.5_dropFalse_rotateFalse_reampleFalse_row1_col0.h5ad')
# # print(sliceC)
#
# pi, log = paste.pairwise_align(sliceA, sliceB, return_obj=True)
#
# print(pi.shape)
# print(np.sum(np.dot(pi, np.ones((pi.shape[1], 1)))))
# print(log)
# print(pi[0][0])


"""
Code for evaluating the scale of gene expression distance
"""

# sliceA_filename = "/Users/xinhaoliu/Desktop/Research/Data/PASTE/delta1/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_delta1_row0_col0.h5ad"
# sliceB_filename = "/Users/xinhaoliu/Desktop/Research/Data/PASTE/delta1/151674_overlap1.5_dropFalse_rotateFalse_reampleTrue_delta1_row0_col1.h5ad"
# dissimilarity = 'kl'
# use_rep = None
#
#
#
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])
#
# spotnamesA = sliceA.obs.index
# spotnamesB = sliceB.obs.index
# common_spots = intersect(spotnamesA, spotnamesB)
# matched_spots = []
# for spot in common_spots:
#     matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))
#
#
# # subset for common genes
# common_genes = intersect(sliceA.var.index, sliceB.var.index)
# sliceA = sliceA[:, common_genes]
# sliceB = sliceB[:, common_genes]
#
# # Calculate expression dissimilarity
# A_X, B_X = to_dense_array(extract_data_matrix(sliceA, use_rep)), to_dense_array(extract_data_matrix(sliceB, use_rep))
# if dissimilarity.lower() == 'euclidean' or dissimilarity.lower() == 'euc':
#     M = distance_matrix(A_X, B_X)
# else:
#     s_A = A_X + 0.01
#     s_B = B_X + 0.01
#     M = kl_divergence(s_A, s_B)
#
# print(M.shape)
# print(len(matched_spots))
# for matched_spot in matched_spots:
#     print("-------------")
#     print(matched_spot)
#     source_spot_idx = matched_spot[0]
#     dest_spot_idx = matched_spot[1]
#     print(M[source_spot_idx][dest_spot_idx])
#     print(M[source_spot_idx].min())

# print(M.shape)
# print(len(matched_spots))
# avg_distance = 0
# for matched_spot in matched_spots:
#     # print("-------------")
#     # print(matched_spot)
#     source_spot_idx = matched_spot[0]
#     dest_spot_idx = matched_spot[1]
#     sum_distance = 0
#     for cost in M[source_spot_idx]:
#         sum_distance += cost - M[source_spot_idx][dest_spot_idx]
#     avg_distance += (sum_distance / len(M[source_spot_idx]))
#
# avg_distance /= len(matched_spots)
# print(avg_distance)

"""
Code for evaluating the scale of gene expression distance ends.
"""


def calculate_overlap_region(sliceA_filename, sliceB_filename):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])

    spotnamesA = sliceA.obs.index
    spotnamesB = sliceB.obs.index
    common_spots = intersect(spotnamesA, spotnamesB)
    matched_spots = []
    for spot in common_spots:
        matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))

    maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
    return maximum_possible_accuracy


sliceA = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151510/151510_overlap0.3_dropFalse_rotateFalse_resampleTrue_delta1.0_row0_col0.h5ad"
sliceB = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151510/151510_overlap0.3_dropFalse_rotateFalse_resampleTrue_delta1.0_row1_col0.h5ad"
overlap = calculate_overlap_region(sliceA, sliceB)
print(overlap)
