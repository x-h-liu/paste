import scanpy as sc
import numpy as np


# adata = sc.read_10x_h5('/n/fs/ragr-data/datasets/RongFan/CUT_RNA_H3K4me3_15/CUT15/raw_peak_bc_matrix.h5')
# print(adata)

# adata = sc.read_10x_h5('/n/fs/ragr-data/datasets/RongFan/CUT_RNA_H3K27ac_14/CUT14/raw_peak_bc_matrix.h5')
# print(adata)

# adata = sc.read_10x_h5('/n/fs/ragr-data/datasets/RongFan/CUT_RNA_H3K27me3_13/CUT13/raw_peak_bc_matrix.h5')
# print(adata)

adata = sc.read_10x_h5('/n/fs/ragr-data/datasets/RongFan/ATAC_RNA/ATAC/raw_peak_bc_matrix.h5')
print(adata)


matrix = np.load("/n/fs/ragr-data/datasets/RongFan/ATAC_RNA/ATAC_lsi.npy")
print(matrix)
print(matrix.shape)
