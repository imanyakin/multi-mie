import sys
from mayavi import mlab

PATH = '/home/ilya/workspace/plasmonics/py_gmm/'
sys.path.append(PATH)
import numpy as np
import numpy.ma as ma
from matplotlib.patches import Circle, Ellipse
import py_gmm
from nplab.modelling.wavelength_to_rgb import wavelength_to_rgb
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from plotly import tools
from mpl_toolkits.mplot3d import Axes3D

from plot_functions import plot_3d_positions,plot_scatter_3d_wireframe,autoscale_equal_axis,plot_scatter_3d_color_sphere
from scattering import 
np.random.seed(3)
DEBUG = 0

'''Globals'''
RADIUS = 30.0
SPACING = 1.0
Ns = 2
eps_db_out=py_gmm.mat.generate_eps_db(PATH+'epsilon/',ext='*.edb')
eps_files,eps_names,eps_db=eps_db_out['eps_files'],eps_db_out['eps_names'],eps_db_out['eps_db']
target_comp= np.array(['eAuJCSTDf']*Ns) # vector containing the optical constants names
n_matrix = 1.33  # water
# Euler angles: (alpha,beta,gamma)=(0,0,0) means a z-directed, x-polarized plane wave


alpha = 0.0   # azimuth
beta = 0.0  # polar
gamma = 0.0  # polarization

print alpha,beta,gamma
d,p =euler_to_plane_wave(alpha,beta,gamma)
print d, p

plane_wave_normal = d 
polarization = p


print "Plane wave propagation direction:", d
print "Polarization", p



if __name__ == "__main__":
	print "pass"