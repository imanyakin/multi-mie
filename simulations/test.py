import sys
PATH = '/home/ilya/workspace/plasmonics/multi-mie/lib/'
sys.path.append(PATH)

import numpy as np

from structures import structure_generator, rotate
from plot_functions import plot_3d_positions
from scattering import multiwavelength_scattering_cone
import utils
from nplab.modelling.wavelength_to_rgb import wavelength_to_rgb

knight_xyz,knight_rs = structure_generator("knight")
def compute_intensity(args):
	(angle, wavelengths) = args
	xyz = rotate(knight_xyz,0,0,angle)
	intensities = multiwavelength_scattering_cone(xyz,knight_rs,wavelengths,debug = 0)
	return intensities
I = []

angles = np.linspace(0,360,36)
wavelengths = [500,550,600,650,700,750,800]

from multiprocessing import Pool

p = Pool()
I = p.map(compute_intensity, [(angle,wavelengths) for angle in angles])


import matplotlib.pyplot as plt 
I = np.array(I)
print I.shape
for i in range(len(wavelengths)):
	wl = wavelengths[i]
	wl_intensity = np.array(I[:,i])
	wl_intensity = wl_intensity / np.max(wl_intensity)

	plt.plot(angles,wl_intensity,"x-",label="Wavelength: {}".format(wl))
plt.legend()
plt.show()
print I.shape