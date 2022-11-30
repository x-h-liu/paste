import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors
import matplotlib.patches as mpatches




"""
slide-seq
"""
# cell_to_color_map = {}
# cell_to_color_map['neuroplacodal cell'] = sns.color_palette('tab20')[0]
# cell_to_color_map['paraxial cell'] = sns.color_palette('tab20')[1]
# cell_to_color_map['native cell'] = sns.color_palette('tab20')[2]
# cell_to_color_map['premigratory neural crest cell'] = sns.color_palette('tab20')[3]
# cell_to_color_map['midbrain dopaminergic neuron'] = sns.color_palette('tab20')[4]
# cell_to_color_map['primitive red blood cell'] = sns.color_palette('tab20')[5]
# cell_to_color_map['mesodermal cell'] = sns.color_palette('tab20')[6]
# cell_to_color_map['neurectodermal cell'] = sns.color_palette('tab20')[7] # E8.5 specific cell type
# cell_to_color_map['endodermal cell'] = sns.color_palette('tab20')[8]
# cell_to_color_map['spinal cord interneuron'] = sns.color_palette('tab20')[9]
# cell_to_color_map['heart valve cell'] = sns.color_palette('tab20')[10]
# cell_to_color_map['surface ectodermal cell'] = sns.color_palette('tab20')[11]
# cell_to_color_map['hemangioblast'] = sns.color_palette('tab20')[12]
# cell_to_color_map['gut endothelial cell'] = sns.color_palette('tab20')[13]
# cell_to_color_map['splanchnic mesodermal cell'] = sns.color_palette('tab20')[14]
# cell_to_color_map['notochordal cell'] = sns.color_palette('tab20')[15] # E9.5 specific cell type
# def plot_slice(slice_,figsize=None,obs_name='cell_type',color_map=cell_to_color_map,ax=None, s=100):
#     (min_x,min_y),(max_x,max_y) = slice_.obsm['X_spatial'].min(axis=0),slice_.obsm['X_spatial'].max(axis=0)
#     len_x,len_y=max_x-min_x,max_y-min_y
#     if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
#     if not ax: plt.figure(figsize=figsize)
#     colors = list(slice_.obs[obs_name].astype('str').map(color_map))
#     g = sns.scatterplot(x = slice_.obsm['X_spatial'][:,0],y = slice_.obsm['X_spatial'][:,1],linewidth=0,s=s, marker=".",c=colors,ax=ax)
#     #plt.legend(handles=[mpatches.Patch(color=cell_to_color_map[slice_.obs[obs_name].cat.categories[i]], label=slice_.obs[obs_name].cat.categories[i]) for i in range(len(slice_.obs[obs_name].cat.categories))],fontsize=10,title='Cell type',title_fontsize=15)
#     if not ax: ax=g
#     if ax:
#         ax.invert_yaxis()
#         ax.axis('off')




# # def plot_slice(slice_,figsize=None,ax=None, s=100):
# #     (min_x,min_y),(max_x,max_y) = slice_.obsm['X_spatial'].min(axis=0),slice_.obsm['X_spatial'].max(axis=0)
# #     len_x,len_y=max_x-min_x,max_y-min_y
# #     if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
# #     if not ax: plt.figure(figsize=figsize)
# #     g = sns.scatterplot(x=slice_.obsm['X_spatial'][:,0],y=slice_.obsm['X_spatial'][:,1],linewidth=0,s=s, marker=".",ax=ax)
# #     if not ax: ax=g
# #     if ax:
# #         ax.invert_yaxis()
# #         ax.axis('off')


# def plot_slice_umi(slice_, figsize=None, ax=None, s=100):
#     (min_x,min_y),(max_x,max_y) = slice_.obsm['X_spatial'].min(axis=0), slice_.obsm['X_spatial'].max(axis=0)
#     len_x, len_y = max_x-min_x, max_y-min_y
#     if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
#     if not ax: plt.figure(figsize=figsize)
#     min_umi = np.min(slice_.obs["nCount_RNA"])
#     max_umi = np.max(slice_.obs["nCount_RNA"])
#     norm = colors.Normalize(vmin=min_umi, vmax=max_umi, clip=True)
#     # cmap = sns.color_palette("rocket_r", as_cmap=True)
#     cmap = sns.color_palette("YlOrBr", as_cmap=True)
#     c = [cmap(norm(i)) for i in slice_.obs["nCount_RNA"]]
#     g = sns.scatterplot(x = slice_.obsm['X_spatial'][:,0],y = slice_.obsm['X_spatial'][:,1],linewidth=0,s=s, marker=".", c=c, ax=ax)
#     if not ax: ax=g
#     if ax:
#         ax.invert_yaxis()
#         ax.axis('off')


# def slice_umi_dist(slice, outlier=0):
#     spot_umi_counts = slice.obs['nCount_RNA']
#     to_keep = (spot_umi_counts <= np.partition(spot_umi_counts, -(outlier + 1))[-(outlier + 1)])
#     slice = slice[slice.obs.index[to_keep]]
#     spot_umi_counts = slice.obs['nCount_RNA']

#     print("mean UMI per spot is: " + str(np.mean(spot_umi_counts)))
#     print("median UMI per spot is: " + str(np.median(spot_umi_counts)))
#     print("max UMI per spot is: " + str(np.max(spot_umi_counts)))
#     print("min UMI per spot is: " + str(np.min(spot_umi_counts)))

#     counts, edges, bars = plt.hist(spot_umi_counts)
#     plt.bar_label(bars)
#     plt.xlabel("UMI")
#     plt.ylabel("Number of spots")



# adata1 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e85_201104_07.h5ad")
# # adata2 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e85_201104_27.h5ad")
# # adata3 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e85_201104_28.h5ad")
# # adata4 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e85_201104_29.h5ad")
# # adata1 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e95_201112_03.h5ad")
# # adata2 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e95_201112_04.h5ad")
# # adata3 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e95_201112_05.h5ad")
# # adata4 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e95_201112_32.h5ad")
# # adata5 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e95_201112_33.h5ad")
# # adata6 = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/slideseq/mouse_embryo/e95_201112_36.h5ad")
# print(adata1)
# # print(adata.obs['cell_type'])
# # print(np.array(adata.obs["cell_type"].unique()))
# # # print(adata.obs['nCount_RNA'])
# # print(adata.obs['nFeature_RNA'])

# slice_umi_dist(adata1)
# plt.show()

# # plot_slice(adata1)
# # plot_slice(adata2)
# # plot_slice(adata3)
# # plot_slice(adata4)
# # plot_slice_umi(adata)
# # plot_slice_umi(adata1)
# # plot_slice_umi(adata2)
# # plot_slice_umi(adata3)
# # plot_slice_umi(adata4)
# # plot_slice_umi(adata5)
# # plot_slice_umi(adata6)
# # plt.show()



"""
stereo-seq
"""
# adata = sc.read_h5ad("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h_a_count_normal_stereoseq.h5ad")
# print(adata)
# print(np.array(adata.obs['annotation'].unique()))


# adata1 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S01']
# adata1.obsm['spatial'] = adata1.obsm['spatial'][:,:2]
# adata1.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s01.h5ad')

# adata2 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S02']
# adata2.obsm['spatial'] = adata2.obsm['spatial'][:,:2]
# adata2.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s02.h5ad')

# adata3 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S03']
# adata3.obsm['spatial'] = adata3.obsm['spatial'][:,:2]
# adata3.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s03.h5ad')

# adata4 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S04']
# adata4.obsm['spatial'] = adata4.obsm['spatial'][:,:2]
# adata4.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s04.h5ad')

# adata5 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S05']
# adata5.obsm['spatial'] = adata5.obsm['spatial'][:,:2]
# adata5.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s05.h5ad')

# adata6 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S06']
# adata6.obsm['spatial'] = adata6.obsm['spatial'][:,:2]
# adata6.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s06.h5ad')

# adata7 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S07']
# adata7.obsm['spatial'] = adata7.obsm['spatial'][:,:2]
# adata7.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s07.h5ad')

# adata8 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S08']
# adata8.obsm['spatial'] = adata8.obsm['spatial'][:,:2]
# adata8.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s08.h5ad')

# adata9 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S09']
# adata9.obsm['spatial'] = adata9.obsm['spatial'][:,:2]
# adata9.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s09.h5ad')

# adata10 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S10']
# adata10.obsm['spatial'] = adata10.obsm['spatial'][:,:2]
# adata10.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s10.h5ad')

# adata11 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S11']
# adata11.obsm['spatial'] = adata11.obsm['spatial'][:,:2]
# adata11.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s11.h5ad')

# adata12 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S12']
# adata12.obsm['spatial'] = adata12.obsm['spatial'][:,:2]
# adata12.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s12.h5ad')

# adata13 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S13']
# adata13.obsm['spatial'] = adata13.obsm['spatial'][:,:2]
# adata13.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s13.h5ad')

# adata14 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S14']
# adata14.obsm['spatial'] = adata14.obsm['spatial'][:,:2]
# adata14.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s14.h5ad')

# adata15 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S15']
# adata15.obsm['spatial'] = adata15.obsm['spatial'][:,:2]
# adata15.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s15.h5ad')

# adata16 = adata[adata.obs['slice_ID'] == 'E14-16h_a_S16']
# adata16.obsm['spatial'] = adata16.obsm['spatial'][:,:2]
# adata16.write('/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s16.h5ad')

def plot_slice_umi(slice_, figsize=None, ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0), slice_.obsm['spatial'].max(axis=0)
    len_x, len_y = max_x-min_x, max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    min_umi = np.min(slice_.obs["sum_umi"])
    max_umi = np.max(slice_.obs["sum_umi"])
    norm = colors.Normalize(vmin=min_umi, vmax=max_umi, clip=True)
    # cmap = sns.color_palette("rocket_r", as_cmap=True)
    cmap = sns.color_palette("YlOrBr", as_cmap=True)
    c = [cmap(norm(i)) for i in slice_.obs["sum_umi"]]
    g = sns.scatterplot(x = slice_.obsm['spatial'][:,0],y = slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".", c=c, ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')


def slice_umi_dist(filename):
    print(filename)
    adata = sc.read_h5ad(filename)
    gene_expression_matrix = adata.layers['raw_counts']
    spot_umi_counts = np.sum(gene_expression_matrix, axis=1)

    print("number of spots is: " + str(len(spot_umi_counts)))
    print("mean UMI per spot is: " + str(np.mean(spot_umi_counts)))
    print("median UMI per spot is: " + str(np.median(spot_umi_counts)))


slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s01.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s02.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s03.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s04.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s05.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s06.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s07.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s08.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s09.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s10.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s11.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s12.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s13.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s14.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s15.h5ad")
slice_umi_dist("/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s16.h5ad")

