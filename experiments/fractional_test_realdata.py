from src.paste.fractional_align import partial_pairwise_align
from experiments.helper import compute_alignment_ari
from experiments.explore_partial_properties import plot_slice_pairwise_alignment
import scanpy as sc
import numpy as np
import pandas as pd
import scipy
import src.paste.PASTE as paste


def mapping_accuracy(labels1, labels2, pi):
    mapping_dict = {'Layer1': 1, 'Layer2': 2, 'Layer3': 3, 'Layer4': 4, 'Layer5': 5, 'Layer6': 6, 'WM': 7}
    return np.sum(pi*(scipy.spatial.distance_matrix(np.matrix(labels1.map(mapping_dict)).T, np.matrix(labels2.map(mapping_dict)).T)==0))


def max_accuracy(labels1, labels2):
    w = min(1/len(labels1), 1/len(labels2))
    cats = set(pd.unique(labels1)).union(set(pd.unique(labels1)))
    return sum([w * min(sum(labels1==c), sum(labels2==c)) for c in cats])

# sliceA_filename = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151508.h5ad"
# sliceB_filename = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151509.h5ad"
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
#
# # pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.05, m=0.99, armijo=False, dissimilarity='kl', norm=True, return_obj=True, verbose=True)
# pi, log = paste.pairwise_align(sliceA, sliceB, alpha=0.05, dissimilarity='kl', norm=True, return_obj=True, verbose=True)
# acc = mapping_accuracy(sliceA.obs['layer_guess_reordered'], sliceB.obs['layer_guess_reordered'], pi)
# print(acc)




# pairs = [(151507, 151508), (151508, 151509), (151509, 151510),
#          (151669, 151670), (151670, 151671), (151671, 151672),
#          (151673, 151674), (151674, 151675), (151675, 151676)]
# # pair_to_overlap = {
# #     (151507, 151508): 0.99, (151508, 151509): 0.5, (151509, 151510): 0.99,
# #     (151669, 151670): 0.99, (151670, 151671): 0.8, (151671, 151672): 0.99,
# #     (151673, 151674): 0.95, (151674, 151675): 0.95, (151675, 151676): 0.99
# # }
# pair_to_overlap = {
#     (151507, 151508): 0.9262388673213062, (151508, 151509): 0.7643125783535312, (151509, 151510): 0.9458838278311742,
#     (151669, 151670): 0.9292984869325996, (151670, 151671): 0.7886119257086998, (151671, 151672): 0.9042033235581622,
#     (151673, 151674): 0.931774415405777, (151674, 151675): 0.9213204951856945, (151675, 151676): 0.9206393718452048
# }
# print("NEW OVERLAP")
# for pair in pairs:
#     print("=========================================")
#     sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/' + str(pair[0]) + '.h5ad'
#     sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/' + str(pair[1]) + '.h5ad'
#     sliceA = sc.read_h5ad(sliceA_filename)
#     sliceB = sc.read_h5ad(sliceB_filename)
#     # print("ORIGINAL")
#     # pi, log = paste.pairwise_align(sliceA, sliceB, alpha=0.9, dissimilarity='kl', norm=True, return_obj=True, verbose=True)
#     print("PARTIAL")
#     pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.9, m=pair_to_overlap[pair], armijo=False, dissimilarity='kl', norm=True,
#                                      return_obj=True, verbose=True)
#     acc = mapping_accuracy(sliceA.obs['layer_guess_reordered'], sliceB.obs['layer_guess_reordered'], pi)
#     max_acc = max_accuracy(sliceA.obs['layer_guess_reordered'], sliceB.obs['layer_guess_reordered'])
#     print(pair)
#     print("Accuracy is " + str(acc))
#     print("Max accuracy is " + str(max_acc))



# sliceA_filename = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151508.h5ad"
# sliceB_filename = "/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151509.h5ad"
sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151670.h5ad'
sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151671.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)

partial_pi, partial_log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=0.5, armijo=False, dissimilarity='glmpca', norm=True, return_obj=True, verbose=True)
full_pi, full_log = paste.pairwise_align(sliceA, sliceB, alpha=0.1, dissimilarity='kl', norm=True, return_obj=True, verbose=True)
partial_ari = compute_alignment_ari(sliceA, sliceB, partial_pi)
full_ari = compute_alignment_ari(sliceA, sliceB, full_pi)
print("partial ari is: " + str(partial_ari))
print("full ari is: " + str(full_ari))
plot_slice_pairwise_alignment(sliceA, sliceB, partial_pi)


