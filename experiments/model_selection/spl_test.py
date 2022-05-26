import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from src.paste.fractional_align import partial_pairwise_align_given_cost_matrix
from src.paste.helper import intersect, to_dense_array, extract_data_matrix, glmpca_distance


def wloss(M, T):
    return np.sum(M * T)


m_to_run = [0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]

# sliceA_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151673/151673_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0.0_row0_col0.h5ad'
# sliceB_filename = '/Users/xinhaoliu/Desktop/Research/Code/st_overlap_sim/sim/151674/151674_overlap0.7_dropFalse_rotateFalse_resampleFalse_delta0.0_row1_col0.h5ad'
# sliceA = sc.read_h5ad(sliceA_filename)
# sliceB = sc.read_h5ad(sliceB_filename)
# common_genes = intersect(sliceA.var.index, sliceB.var.index)
# sliceA = sliceA[:, common_genes]
# sliceB = sliceB[:, common_genes]
# # Get transport cost matrix
# A_X, B_X = to_dense_array(extract_data_matrix(sliceA, None)), to_dense_array(extract_data_matrix(sliceB, None))
# M = glmpca_distance(A_X, B_X, latent_dim=50, filter=True)
#
# w_losses = []
# for m in m_to_run:
#     print("============================")
#     print("m = " + str(m))
#     pi, log = partial_pairwise_align_given_cost_matrix(sliceA, sliceB, M=M, alpha=0.1, m=m, armijo=False,
#                                                        norm=True, return_obj=True, verbose=True)
#     w_losses.append(wloss(M, pi))
#
# print(w_losses)


def plot_cost_curve(m_list, cost_list):
    fig, ax = plt.subplots()
    ax.plot(m_list, cost_list)
    ax.set_xlim(1, 0)
    ax.set_xticks([0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])
    ax.set_xlabel("m")
    ax.set_ylabel("W loss")
    # ax.set_title("Delta=1.0, Truth=0.9")
    ax.set_title("real0010, Truth=0.7")
    plt.show()


w_losses = [20.77169504621095, 19.062462147642695, 17.188081809248306, 15.600164363664597, 14.129847262734824, 12.839036424127363, 11.65402368370547, 10.54504705444226, 9.50615487289093, 8.522570194847015, 7.577489260327491, 6.669950795207466, 5.793752240175898, 4.947794128975601, 4.1410079454594815, 3.3601879158865215, 2.611958737181231, 1.895583017548707, 1.211870604510512, 0.5694889094752186]
plot_cost_curve(m_to_run, 100 * np.array(w_losses))





# pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=0.8, armijo=False, dissimilarity='glmpca', norm=True, return_obj=True, verbose=True)
# print(pi.shape)
#
# going_out_mass_A = np.sum(pi, axis=1)
# total_mass_A = np.ones((sliceA.shape[0],)) / sliceA.shape[0]
# sunk_mass_A = total_mass_A - going_out_mass_A
# # average_mass_A = going_out_mass_A[np.nonzero(going_out_mass_A)].mean()
# # max_mass_A = 1.0 / sliceA.shape[0]
# # print("average mass going out per spot: " + str(average_mass_A))
# # print("max possible mass going out per spot: " + str(max_mass_A))
#
#
#
# counts, edges, bars = plt.hist(sunk_mass_A, bins=20)
# plt.bar_label(bars)
# plt.xlabel("Sunk mass")
# plt.ylabel("Number of spots")
# plt.show()
