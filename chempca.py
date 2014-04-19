import numpy as np
from scipy import linalg as LA
import astropy.io.fits as fits
import matplotlib.pyplot as p

def PCA(data, dims_rescaled_data=2):
    """
    returns: data transformed in 2 dims/columns + regenerated original data
    pass in: data as 2D NumPy array
    """
    mn = np.mean(data, axis=0)
    # mean center the data
    data -= mn
    # calculate the covariance matrix
    C = np.cov(data.T)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    evals, evecs = LA.eig(C)
    # sorted them by eigenvalue in decreasing order
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:,idx]
    evals = evals[idx]
    print(evals)
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or dims_rescaled_data)
#    evecs = evecs[:,:dims_rescaled_data]
    # carry out the transformation on the data using eigenvectors
    data_rescaled = np.dot(evecs.T, data.T).T
    # reconstruct original data array
#    data_original_regen = np.dot(evecs,data).T + mn
    return data_rescaled,evals#, data_original_regen


cprops = fits.getdata('multitracer.fits')


array = np.log10(np.vstack((cprops['HCN'],cprops['HCOP'],cprops['H13COP'],cprops['H13CN'],cprops['HNCO'])))#,cprops['HNCO'],cprops['CS'],cprops['SIO'])))
array[np.isnan(array)]=np.nanmin(array)

rs,evals = PCA(array.T)



