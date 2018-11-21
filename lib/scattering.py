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


def euler_to_plane_wave(alpha,beta,gamma,units="degrees"):
	if units == "degrees":
		alpha = alpha*(np.pi/180.0)
		beta = beta*(np.pi/180.0)
		gamma = gamma*(np.pi/180.0)
	dx = np.cos(alpha)*np.sin(beta)
	dy = np.sin(alpha)*np.sin(beta)
	dz = np.cos(beta)

	Ex = np.cos(alpha)*np.cos(beta)*np.cos(gamma)-np.sin(alpha)*np.sin(gamma)
	Ey = np.sin(alpha)*np.cos(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma)
	Ez =-np.sin(beta)*np.cos(gamma)
	if np.absolute(Ez) < 1e-12:
		Ez = 0.0 
	assert(np.dot([dx,dy,dz],[Ex,Ey,Ez])==0)
	return [dx,dy,dz], [Ex,Ey,Ez]

def define_incident(debug = 0):

	alpha = 0.0
	beta = 0.0
	gamma = 0.0

	if debug > 0:
		direction, E = euler_to_plane_wave(alpha,beta,gamma)
		print "Light propagation direction", direction
		print "Light polarization", E
		#test that polarization is perp to propagation
		assert(np.dot(direction,E) == 0)

	return alpha,beta,gamma

def simulation_configuration(n_particle):
	eps_db_out=py_gmm.mat.generate_eps_db(PATH+'epsilon/',ext='*.edb')
	eps_files,eps_names,eps_db=eps_db_out['eps_files'],eps_db_out['eps_names'],eps_db_out['eps_db']
	target_comp= np.array(['eAuJCSTDf']*n_particle) # vector containing the optical constants names
	n_matrix = 1.33  # water

	return eps_db, target_comp, n_matrix

def simulate_scattering(particle_xyz,wavelength,particle_radii):
	
	''' Get configuration for simulation '''
	n_particle = particle_xyz.shape[0]
	eps_db, target_comp, n_matrix = simulation_configuration(n_particle)
	alpha,beta, gamma = define_incident()


	''' Simulate configuration '''

	n_stop=10 # maximum multipolar expansion order
	f_int=0.0; # interaction cutoff (normally set to zero to have full interactions)
	qs_flag='no' # retarded simulation
	e_list=py_gmm.mat.db_to_eps(wavelength,eps_db,target_comp)
	m_eps=np.column_stack((np.real(e_list),np.imag(e_list)))
	out=py_gmm.gmm_py.gmm_f2py_module.expansion_coefficients(particle_xyz,  # target sphere position in nm
		particle_radii,  # target sphere radii in nm
		m_eps,  # e1 and e2 for each sphere
		f_int,  # interaction coefficient
		n_matrix,  # environment refractive index
		wavelength, # computation wavelength
		alpha,beta,gamma,  # euler angles for the incident pw
		n_stop,  # maximum number for expansion coefficients
		qs_flag)  # quasi static approximation
	v_amnbmn=out[0][:,0] # getting field expansion coefficients
	v_dmncmn=out[0][:,1]

    # local field
	v_emn=py_gmm.gmm_py.gmm_f2py_module.emn(n_stop)[0] # normalization coeffs

    # retrieving the scattered cloud
	k_far = 2.0*np.pi*n_matrix/wavelength
	m_sc, scatot, error = py_gmm.gmm_py.gmm_f2py_module.efar(k_far,
		n_stop, #ns
		0.0,360.0,360, #phimin,phimax,phistep
		180, #theta_step
		0.0, #betap
		v_amnbmn, #v_amnbmn
		particle_xyz) 
	return m_sc

def multiwavelength_scattering(particle_xyz,particle_radii,wavelengths,debug=1):
	if debug > 0:
		ax=plot_3d_positions(particle_xyz,particle_radii,with_mayavi=True)
	for wl in wavelengths:
		m_sc_out = simulate_scattering(particle_xyz,wl,particle_radii)
		if debug > 0:
			color = [c/255.0 for c in wavelength_to_rgb(wl)]
			plot_scatter_3d_wireframe(m_sc_out,ax=ax,color=color,label="Wavelength:{0} nm".format(wl),with_mayavi=True)
	# autoscale_equal_axis(ax)
	# ax.legend()

def multiwavelength_scattering_difference(particle_xyz,particle_radii,wl_1,wl_2,debug=1):
	if debug > 0:
		ax=plot_3d_positions(particle_xyz,particle_radii)

	m_sc_out_1 = simulate_scattering(particle_xyz,wl_1,particle_radii)
	m_sc_out_2 = simulate_scattering(particle_xyz,wl_2,particle_radii)

	theta = m_sc_out_1[:,:,1].real
	phi = m_sc_out_1[:,:,0].real
	r_1 = m_sc_out_1[:,:,6].real
	r_2 = m_sc_out_2[:,:,6].real
    
	scalefactor = min(np.max(r_1),np.max(r_2))
	r_1 = r_1 / np.max(r_1)
	r_2 = r_2 / np.max(r_2)
    
	r_diff = np.absolute(r_1 - r_2)

	x = scalefactor*r_diff * np.sin(phi) * np.cos(theta)
	y = scalefactor*r_diff * np.cos(phi)
	z = scalefactor*r_diff * np.sin(phi) * np.sin(theta)

	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')
	ax.plot_wireframe(x,y,z,alpha=0.05)
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")

def multiwavelength_scattering_cone(particle_xyz,particle_radii,wavelengths,ax=None,debug=1):

	if ax is None and debug > 0:
		fig = plt.figure()
		ax=plot_3d_positions(particle_xyz,particle_radii)
	intensities = []

	for wl in wavelengths:
		m_sc_out = simulate_scattering(particle_xyz,wl,particle_radii)
		theta = m_sc_out[:,:,1].real #in xy plane
		phi = m_sc_out[:,:,0].real #from z-axis
		r = m_sc_out[:,:,6].real

		phi_down = 180*np.pi/180.0
		phi_center = 180*np.pi/180.0
		phi_up = 0*np.pi/180.0
		phi_bounds = [max(0.0,phi_center-phi_down),min(np.pi,phi_center+phi_up)]
		phi_inds = np.logical_and(phi >= phi_bounds[0], phi <= phi_bounds[1])

		theta_lower = 0.0
		theta_upper = 2*np.pi
		theta_bounds= [theta_lower,theta_upper]
		theta_inds = np.logical_and(theta >= theta_bounds[0], theta <= theta_bounds[1])
		
		inds = np.logical_and(phi_inds,theta_inds)

	
		r_new = np.zeros(theta.shape)
		r_new[inds] = r[inds]
		r = r_new

		#phi - inclanation relative to z-axis
		#theta - scattering angle relative to x-axis?
		x = r * np.sin(phi) * np.cos(theta)
		y = r * np.sin(phi) * np.sin(theta)
		z = r * np.cos(phi)

		total = np.sum(np.sqrt(x**2 + y**2 + z**2))

		if debug > 0:
			color = [c/255.0 for c in wavelength_to_rgb(wl)]
			ax.plot_wireframe(x,y,z,alpha=0.5,color=color,label="Wavelength {} nm".format(wl))
			ax.set_xlim(-2000,2000)
			ax.set_ylim(-2000,2000)
			ax.set_zlim(-2000,2000)

		intensities.append(total)
	return intensities