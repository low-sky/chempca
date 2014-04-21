import astropy.io.fits
from spectral_cube.spectral_cube import SpectralCube, SpectralCubeMask, LazyMask, FunctionMask, read
from signal_id import noise
input_file = '/Users/erosolo/code/noodle/n253_hcn.fits'
cube = read(input_file,format='fits')
noiseobj = noise.Noise(cube)
noiseobj.mask_out_signal()
