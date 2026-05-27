import numpy as np

X_orig = np.loadtxt('X.csv', dtype = np.float32, skiprows = 1, delimiter = ',')


# standardize the data
X_centered = (X_orig - np.mean(X_orig, axis = 0)) / np.std(X_orig, axis = 0)
X_centered = X_centered.astype(np.float64)


# doing PCA
cov_mat = np.cov(X_centered, rowvar = False)
eigenvalues, eigenvectors = np.linalg.eigh(cov_mat)


# sorting the eigenvalues and eigenvectors
sort_idx = np.argsort(eigenvalues)[::-1]
eigenvalues_sorted = eigenvalues[sort_idx]
eigenvectors_sorted = eigenvectors[:, sort_idx]

explained_variance = eigenvalues_sorted/np.sum(eigenvalues_sorted)
cumulative_explained_variance = np.cumsum(explained_variance)

k = 1
for i in range(0, cumulative_explained_variance.shape[0]):
    
    prop_expt_var = cumulative_explained_variance[i]
    
    if prop_expt_var >= 0.990:  # select the top k components that explain at least 99% of the total variance
        break

    k = k + 1



print(f'selecting {k} components')


# do the dimensionality reduction
X_reduced = X_centered @ eigenvectors_sorted[:,:k]

np.savetxt('X_pca_reduced.csv', X_reduced, delimiter = ',', fmt = '%.4f')

print(X_orig.shape)
print(X_reduced.shape)