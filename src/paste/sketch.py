import numpy
import numpy as np
from scipy.stats import entropy
import scanpy as sc
# from src.paste.helper import to_dense_array, extract_data_matrix
# from helper import to_dense_array, extract_data_matrix
# from src.paste.fractional_align import partial_pairwise_align
# from experiments.helper import plot_slice_pairwise_alignment


# gene_umi_counts = np.array([0, 5, 4, 3, 3, 2, 7, 9])
# print(np.sort((-gene_umi_counts).argsort()[:3]))


# # sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674.h5ad'
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/151674.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# gene_expression_matrix = to_dense_array(extract_data_matrix(sliceA, None))
#
# print(gene_expression_matrix.shape)
# spot_umi_counts = np.sum(gene_expression_matrix, axis=1)
# print(spot_umi_counts.shape)
# print("mean UMI per spot is: " + str(np.mean(spot_umi_counts)))
# print("median UMI per spot is: " + str(np.median(spot_umi_counts)))
# print("max UMI per spot is: " + str(np.max(spot_umi_counts)))
# print("min UMI per spot is: " + str(np.min(spot_umi_counts)))



# sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151670.h5ad'
# sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151671.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
#
# pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=0.7, armijo=False, dissimilarity='glmpca', norm=True, return_obj=True, verbose=True)
#
# plot_slice_pairwise_alignment(sliceA, sliceB, pi)


# probs = [0, 12, 12]
# # h = 0
# # for i in probs:
# #     h += -i * np.log(i)
# #
# # print(h)
# print(entropy(probs))

# lst = ["a", "b", "c", "d"]
# for i in range(len(lst) - 1, -1, -1):
#     print(i)
#     print(lst[i])


# A = np.asarray([[1,2,3],[4,5,6], [7,8,9]])
# # print(np.mean(A, axis=0))
# # print((A - np.mean(A, axis=0)) / np.std(A, axis=0))

# print(A)
# print(A[0,:])
# A += np.array([1, 1, 1])
# print(A)



lst = ['a', 'b', 'c', 'd', 'e']

for i in range(len(lst) - 1, -1, -1):
    print(lst[i])


