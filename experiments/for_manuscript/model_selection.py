import sys
sys.path.insert(0, '/n/fs/ragr-research/users/xinhao/workspace/code/paste')
import scanpy as sc

from src.paste.fractional_model_selection import decide_overlap


"""
Code for testing model selection on simulated data
"""
# overlap_to_run = [0.9, 0.7, 0.5, 0.3]
# delta_to_run = [0.1, 1.0, 2.0, 3.0]

# overlap_delta_estimation = {}
# for overlap in overlap_to_run:
#     overlap_delta_estimation[overlap] = {}
#     for delta in delta_to_run:
#         sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/single_resample/151674_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row0_col0.h5ad'
#         sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/single_resample/151674_overlap' + str(overlap) + '_dropFalse_rotateFalse_resampleTrue_delta' + str(delta) + '_row1_col0.h5ad'
#         sliceA = sc.read_h5ad(sliceA_filename)
#         sliceB = sc.read_h5ad(sliceB_filename)
#         print("=======================")
#         print("Overlap={0}, Delta={1}".format(overlap, delta))
#         estimation = decide_overlap(sliceA, sliceB)
#         print("Estimation of overlap is: " + str(estimation))
#         overlap_delta_estimation[overlap][delta] = estimation


# print("\n****************************************************************************")
# for overlap in overlap_to_run:
#     for delta in delta_to_run:
#         print(overlap, delta)
#         print(overlap_delta_estimation[overlap][delta])



"""
Code for testing model selection on DLPFC top down
"""
print("Patient 1, AB")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151507/twopieces/151507_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151508/twopieces/151508_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))

print("Patient 1, BC")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151508/twopieces/151508_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151509/twopieces/151509_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))

print("Patient 1, CD")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151509/twopieces/151509_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151510/twopieces/151510_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))


print("Patient 2, AB")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151669/twopieces/151669_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151670/twopieces/151670_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))

print("Patient 2, BC")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151670/twopieces/151670_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151671/twopieces/151671_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))

print("Patient 2, CD")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151671/twopieces/151671_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151672/twopieces/151672_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))


print("Patient 3, AB")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151673/twopieces/151673_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151674/twopieces/151674_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))

print("Patient 3, BC")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151674/twopieces/151674_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151675/twopieces/151675_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))

print("Patient 3, CD")
sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151675/twopieces/151675_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151676/twopieces/151676_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row1_col0.h5ad'
sliceA = sc.read_h5ad(sliceA_filename)
sliceB = sc.read_h5ad(sliceB_filename)
print(decide_overlap(sliceA, sliceB))





"""
Code for testing model selection on DLPFC left right
"""
# print("Patient 1, AB")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151507/twopieces_lr/151507_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151508/twopieces_lr/151508_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))

# print("Patient 1, BC")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151508/twopieces_lr/151508_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151509/twopieces_lr/151509_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))

# print("Patient 1, CD")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151509/twopieces_lr/151509_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151510/twopieces_lr/151510_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))


# print("Patient 2, AB")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151669/twopieces_lr/151669_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151670/twopieces_lr/151670_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))

# print("Patient 2, BC")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151670/twopieces_lr/151670_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151671/twopieces_lr/151671_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))

# print("Patient 2, CD")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151671/twopieces_lr/151671_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151672/twopieces_lr/151672_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))


# print("Patient 3, AB")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151673/twopieces_lr/151673_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151674/twopieces_lr/151674_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))

# print("Patient 3, BC")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151674/twopieces_lr/151674_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151675/twopieces_lr/151675_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))

# print("Patient 3, CD")
# sliceA_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151675/twopieces_lr/151675_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col0.h5ad'
# sliceB_filename = '/n/fs/ragr-data/users/xinhao/DLPFC/sim/151676/twopieces_lr/151676_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0_row0_col1.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# print(decide_overlap(sliceA, sliceB))