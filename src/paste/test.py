import numpy as np
import scanpy as sc
import anndata as ad
import src.paste.PASTE as paste


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



p = np.array([[1, 2, 3], [3, 2, 1]])
print(np.argmin(p[1]))


