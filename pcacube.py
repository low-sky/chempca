import numpy as np
#from scipy import linalg as LA
from scipy import stats as ss
import astropy.io.fits as fits
import matplotlib.pyplot as plt

TracerArray = np.array(['cs','h13cn','h13co+',\
                           'hcn','hco+','hnco','hco',\
                           'ch3ocho','ch3oh_extended',\
                           'ch3sh_extended','hn13c','sio','hc3n_extended'])

TracerArray = np.array(['hcn','hco+','hnco','sio'])    
dir = '/srv/astro/erosolo/n253/cubes/newrelease/lines/robust/non_pbcor/'


tracer = 'hcn'
cube= fits.getdata(dir+'ngc253_'+tracer+'_clean_RO.fits')
cube.shape = cube.shape[1:] # Drop polarization axis
hdr = fits.getheader(dir+'ngc253_'+tracer+'_clean_RO.fits')
SignalIndex = np.where(cube > 0.01)

normalization = np.zeros(len(TracerArray))

matrix = None
for idx, tracer in enumerate(TracerArray):
    cube= fits.getdata(dir+'ngc253_'+tracer+'_clean_RO.fits')
    cube.shape = cube.shape[1:]
    vector = cube[SignalIndex]
    vector = ss.rankdata(vector)/len(vector)
    #    score90 = ss.scoreatpercentile(vector,90.0)
    #    vector[np.isnan(vector)]=0.0
    #    normalization[idx] = score90
    #    vector = vector / score90
    if matrix is None:
        matrix = vector
    matrix = np.vstack((matrix,vector))        


        
matrix = matrix.T
means = np.mean(matrix,axis=0)
matrix -= means

Covariance = np.cov(matrix.T)
Evals,Evecs = np.linalg.eig(Covariance)
idx = np.argsort(Evals)[::-1]
Evals = Evals[idx]
Evecs = Evecs[:,idx]
# pull out Evec for max Eval

for Comp,Value in enumerate(Evals):
    wts = Evecs[:,Comp]*np.linalg.det(Evecs)
    WeightCube=np.zeros(cube.shape)
    for index, tracer in enumerate(TracerArray):
        cube= fits.getdata(dir+'ngc253_'+tracer+'_clean_RO.fits')
        cube.shape = cube.shape[1:]
        cube[np.isnan(cube)]=0.0
        WeightCube = wts[index]*cube+WeightCube

    WeightCube[WeightCube==0]=np.nan
    fits.writeto('ngc253.pca_component_'+str(Comp)+'.fits',WeightCube,hdr,clobber=True)



