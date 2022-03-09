import numpy as np
import scanpy as sc
import anndata as ad
import src.paste.PASTE as paste
from src.paste.helper import intersect, to_dense_array, extract_data_matrix, glmpca_distance, pca_distance
import src.paste.helper as helper
from scipy.spatial import distance_matrix
import scipy


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


"""
Code for calculating percentage of overlap between two sub-slices
"""
# def calculate_overlap_region(sliceA_filename, sliceB_filename):
#     sliceA = sc.read_h5ad(sliceA_filename)
#     sliceB = sc.read_h5ad(sliceB_filename)
#     maximum_num_spots = max(sliceA.shape[0], sliceB.shape[0])
#
#     spotnamesA = sliceA.obs.index
#     spotnamesB = sliceB.obs.index
#     common_spots = intersect(spotnamesA, spotnamesB)
#     matched_spots = []
#     for spot in common_spots:
#         matched_spots.append((spotnamesA.get_loc(spot), spotnamesB.get_loc(spot)))
#
#     maximum_possible_accuracy = len(common_spots) / float(maximum_num_spots)
#     return maximum_possible_accuracy
#
#
# sliceA = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151510/151510_overlap0.3_dropFalse_rotateFalse_resampleTrue_delta1.0_row0_col0.h5ad"
# sliceB = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151510/151510_overlap0.3_dropFalse_rotateFalse_resampleTrue_delta1.0_row1_col0.h5ad"
# overlap = calculate_overlap_region(sliceA, sliceB)
# print(overlap)
"""
Code for calculating percentage of overlap between two sub-slices ends
"""


"""
Code for testing KL divergence
"""
# def kl_divergence(X, Y):
#     """
#     Returns pairwise KL divergence (over all pairs of samples) of two matrices X and Y.
#
#     param: X - np array with dim (n_samples by n_features)
#     param: Y - np array with dim (m_samples by n_features)
#
#     return: D - np array with dim (n_samples by m_samples). Pairwise KL divergence matrix.
#     """
#     assert X.shape[1] == Y.shape[1], "X and Y do not have the same number of features."
#
#     X = X / X.sum(axis=1, keepdims=True)
#     Y = Y / Y.sum(axis=1, keepdims=True)
#     log_X = np.log(X)
#     log_Y = np.log(Y)
#     X_log_X = np.matrix([np.dot(X[i], log_X[i].T) for i in range(X.shape[0])])
#     # print(X_log_X.T.shape)
#     # print(np.dot(X, log_Y.T).shape)
#     D = X_log_X.T - np.dot(X, log_Y.T)
#     return np.asarray(D)
#
# def generalized_kl_divergence(X, Y):
#     """
#     Returns pairwise generalized KL divergence (over all pairs of samples) of two matrices X and Y.
#
#     param: X - np array with dim (n_samples by n_features)
#     param: Y - np array with dim (m_samples by n_features)
#
#     return: D - np array with dim (n_samples by m_samples). Pairwise KL divergence matrix.
#     """
#     assert X.shape[1] == Y.shape[1], "X and Y do not have the same number of features."
#     log_X = np.log(X)
#     log_Y = np.log(Y)
#     X_log_X = np.matrix([np.dot(X[i], log_X[i].T) for i in range(X.shape[0])])
#     D = X_log_X.T - np.dot(X, log_Y.T)
#     sum_X = np.sum(X, axis=1)
#     sum_Y = np.sum(Y, axis=1)
#     D = (D.T - sum_X).T + sum_Y.T
#     return np.asarray(D)
#
# def pairwise_align(sliceA, sliceB, use_rep=None):
#
#     # subset for common genes
#     common_genes = intersect(sliceA.var.index, sliceB.var.index)
#     sliceA = sliceA[:, common_genes]
#     sliceB = sliceB[:, common_genes]
#
#     # Calculate expression dissimilarity
#     A_X, B_X = to_dense_array(extract_data_matrix(sliceA, use_rep)), to_dense_array(extract_data_matrix(sliceB, use_rep))
#     s_A = A_X + 0.01
#     s_B = B_X + 0.01
#     M = generalized_kl_divergence(s_A, s_B)
#     M /= M[M > 0].max()
#     M *= 10
#     return M
#     #return kl_divergence(s_A, s_B)
#
# sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row0_col0.h5ad'
# sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/single_resample/151674_overlap0.5_dropFalse_rotateFalse_resampleTrue_delta1.0_row1_col0.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
#
# M = pairwise_align(sliceA, sliceB)
# print(M)
# print(M.max())

"""
Code for testing KL divergence ends
"""



"""
Code for calculating binomial deviance
"""
# def binomial_deviance(x, y):
#     n = np.sum(x)
#     n_prime = np.sum(y)
#
#     # x = x + n / 100.0
#     # y = y + n_prime / 100.0
#     x = x + 0.01
#     y = y + 0.01
#     print(x)
#     print(y)
#
#     deviance = 0
#     for j in range(len(x)):
#         deviance += x[j] * np.log((n_prime * x[j]) / (n * y[j])) + (n - x[j]) * np.log((n - x[j]) / (n * (1 - y[j] / n_prime)))
#     return deviance
#
#
# smallcount_orig = np.array([0, 0, 1, 1, 0, 1])
# smallcount_resample = np.array([1.0/3, 1.0/3, 2.0/3, 2.0/3, 1.0/3, 2.0/3])
#
# largecount_orig = np.array([0, 0, 100, 100, 0, 100])
# largecount_resample = np.array([(1.0/306)*300, (1.0/306)*300, (101.0/306)*300, (101.0/306)*300, (1.0/306)*300, (101.0/306)*300])

# print(binomial_deviance(smallcount_orig, smallcount_resample))
# print(binomial_deviance(largecount_orig, largecount_resample))

"""
Code for calculating binomial deviance ends
"""




"""
Code for testing glmpca
"""
sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151673/151673_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0.0_row0_col0.h5ad'
sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674/151674_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0.0_row0_col1.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
common_genes = intersect(sliceA.var.index, sliceB.var.index)
sliceA = sliceA[:, common_genes]
sliceB = sliceB[:, common_genes]
A_X, B_X = to_dense_array(extract_data_matrix(sliceA, None)), to_dense_array(extract_data_matrix(sliceB, None))


glmpca_distance(A_X, B_X, 20)
# pca_distance(sliceA, sliceB, 5000, 500)
# print(helper.high_umi_gene_distance(A_X, B_X, 2000).shape)
