import numpy as np


X_orig = np.loadtxt('X.csv', dtype = np.float32, skiprows = 1, delimiter = ',')


# stnadaridze the data
X_centered = (X_orig - np.mean(X_orig, axis = 0)) / np.std(X_orig, axis = 0)

# perform SVD
U, S, Vt = np.linalg.svd(X_centered, full_matrices = False)



'''
================================ check the explained variance ===============================================
'''

explained_variance = S**2
total_variance = explained_variance.sum()
explained_variance_ratio = explained_variance/total_variance
cumulative_explained_variance = np.cumsum(explained_variance_ratio)



'''
================== based on the explained variance select the number of components to keep =========================
'''
k = 1
for i in range(0, cumulative_explained_variance.shape[0]):
    
    prop_expt_var = cumulative_explained_variance[i]
    
    if prop_expt_var >= 0.990:  # select the top k components that explain at least 99% of the total variance
        break

    k = k + 1

# for the initial dataset it is 83

print(f'selecting {k} components')


X_reduced = U[:,:k] @ np.diag(S[:k])

np.savetxt('X_svd_reduced.csv', X_reduced, delimiter = ',', fmt = '%.4f')