import numpy
import numpy as np



gene_umi_counts = np.array([0, 5, 4, 3, 3, 2, 7, 9])
print(np.sort((-gene_umi_counts).argsort()[:3]))
