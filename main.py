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

from utils import plot_3d_positions,plot_scatter_3d_wireframe,autoscale_equal_axis,plot_scatter_3d_color_sphere
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


# def plane_wave_to_euler(direction_xyz=[0,0,1],polarization_xyz=[0,1,0]):
# 	'''
# 		v_kinc(1)=k*COS(alpha0)*SIN(beta0)
# 	    v_kinc(2)=k*SIN(alpha0)*SIN(beta0)
# 	    v_kinc(3)=k*COS(beta0)

# 	    !Calcolo E
# 	    v_Einc(1)=COS(alpha0)*COS(beta0)*COS(gamma0)-SIN(alpha0)*SIN(gamma0)
# 	    v_Einc(2)=SIN(alpha0)*COS(beta0)*COS(gamma0)-COS(alpha0)*SIN(gamma0)
# 	    v_Einc(3)=-SIN(beta0)*COS(gamma0)
# 	'''
# 	dx,dy,dz = direction_xyz[0],direction_xyz[1],direction_xyz[2]
# 	Ex,Ey,Ez = polarization_xyz[0],polarization_xyz[1],polarization_xyz[2]
# 	assert(np.dot([dx,dy,dz],[Ex,Ey,Ez])==0)

# 	k = np.sqrt(dx**2 + dy**2 + dz**2)
	
# 	alpha = np.arctan2(dy,dx) #correct arc-tanning, takes account of quadrant
# 	beta = np.arccos(dz)
# 	gamma = np.arctan2((Ex/Ez + np.cos(alpha)/np.tan(beta)),np.sin(alpha))
# 	alpha = alpha*(180.0/np.pi)
# 	beta  = beta *(180.0/np.pi)
# 	gamma = gamma*(180.0/np.pi)
# 	return alpha,beta,gamma

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

def spherical_polar_to_xyz(r,theta,phi,units="radians"):
	'''
	Mathematics convention for angle conversion (note: not ISO convention):
	theta - azimuthal angle
	phi - polar angle 
	r - radius

	'''
	if units == "radians":
		pass		
	elif units == "degrees":
		theta = theta*(180.0/np.pi)
		phi = phi*(180.0/np.pi)
	x = r*np.sin(phi)*np.cos(theta)
	y = r*np.sin(phi)*np.sin(theta)
	z = r*np.cos(phi)

	return {"x":x,"y":y,"z":z}

def xyz_to_spherical_polar(x,y,z):
	'''
	Mathematics convention for angle conversion (note: not ISO convention):
	theta - azimuthal angle
	phi - polar angle 
	r - radius
	'''

	r = np.sqrt(x**2+y**2+z**2)
	phi = np.arccos(z/float(r))
	theta =np.arctan2(y,float(x))

	return {"r":r,"theta":theta,"phi":phi}


def verify_distances(positions):
	dist = RADIUS*100.0
	for i in range(len(positions)-1):
		for j in range(i+1,len(positions)-1):
			p0 = positions[i]
			p1 = positions[j]
			dist = min(dist,np.linalg.norm(p0-p1))
	return dist #minimum distance between two NPs - should be 2*radius

def revalidate(old_positions, new_position):
	valid = False
	current_position = new_position
	while valid == False:
		valid = True
		for p in old_positions:
			r = current_position - p 
			if np.linalg.norm(r) < 2*RADIUS:
				valid = False
				current_position = current_position + r

			if np.linalg.norm(r) == 2*RADIUS:
				valid = False 
				current_position = current_position + r + np.random.uniform(0,1,size=3)
	return current_position

def get_COM(particle_xyz):
	# print "Particle,SHAPE:",particle_xyz.shape
	COM = np.sum(particle_xyz,axis=0)/float(particle_xyz.shape[0])
	# print "COM,SHAPE:",COM.shape
	return COM

def center_COM(particle_xyz):
	com = get_COM(particle_xyz)
	R = np.asarray(particle_xyz.shape[0]*[com])
	particle_xyz = particle_xyz - R
	particle_xyz.shape
	return particle_xyz

def rotate(particles, theta_x,theta_y, theta_z):
	#specified in degrees
	tx = theta_x*np.pi/180.0
	ty = theta_y*np.pi/180.0
	tz = theta_z*np.pi/180.0
	
	# print "tx,ty,tz",tx,ty,tz
	def rotate_x(ps):
		Rx = np.asarray([
			[1.0,     0.0,    0.0          ],
			[0.0, np.cos(tx), -np.sin(tx)  ],
			[0.0, np.sin(tx),  np.cos(tx)  ]
		])

		outp = np.dot(Rx,ps.T).T
		return outp

	def rotate_y(ps):
		Ry = np.asarray([
			[np.cos(ty),  0.0,  np.sin(ty) ],
			[0.0,         1.0,   0.0       ],
			[-np.sin(ty), 0.0,  np.cos(ty) ]
		])

		outp = np.dot(Ry,ps.T).T
		return outp

	def rotate_z(ps):
		Rz = np.asarray([
			[np.cos(tz), -np.sin(tz), 0.0 ],
			[np.sin(tz),  np.cos(tz), 0.0 ],
			[0.0,              0.0,   1.0 ]
		])
		outp = np.dot(Rz,ps.T).T
		return outp

	# print "start:",particles
	outp = rotate_x(particles)
	# print "Rx:",outp
	outp = rotate_y(outp)
	# print "Ry:",outp
	outp = rotate_z(outp)
	# print "Rz:",outp
	return outp

def make_structure(case="linear_chain"):
	particle_radii = np.array([RADIUS]*Ns)
	
	#Generate linear chain
	if case == "linear_chain":
		particle_xyz = np.asarray([[i*(2*RADIUS+SPACING),0.0,0.0] for i in range(Ns) ])
		particle_xyz = center_COM(particle_xyz)
	#Generate random structure
	else:
		particle_xyz = []
		for n in range(Ns):
			if n == 0:
				particle_xyz = particle_xyz + [[0.0,0.0,0.0]]
			else:
				direction = np.random.normal(0,1,size=3)
				direction = direction/np.linalg.norm(direction)
				next_pos = (2*RADIUS+0.9)*direction+np.asarray(particle_xyz[-1])
				particle_xyz = particle_xyz + [revalidate(particle_xyz,next_pos)]
		particle_xyz = np.asarray(particle_xyz)
		min_dist = verify_distances(particle_xyz)
		particle_xyz = center_COM(particle_xyz)
	return particle_xyz,particle_radii

def simulate_scattering(particle_xyz,wavelength,particle_radii):
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

# def cluster_rotations(orig_particle_xyz, angles, wavelength,particle_radii,debug=1):
# 	if debug > 0:
# 		plot_3d_positions(orig_particle_xyz,particle_radii)
# 	intensities = []
# 	for angle in angles:
# 		rotated_particle_xyz = rotate(orig_particle_xyz,0.0,angle,0.0)
# 		if debug > 0:
# 			ax = plot_3d_positions(rotated_particle_xyz,particle_radii)
# 		m_sc_out = simulate_scattering(rotated_particle_xyz,wavelength,particle_radii)
# 		if debug > 0:
# 			plot_scatter_3d_wireframe(m_sc_out,ax=ax,show_plot=True)
# 		v_theta = m_sc_out[:,0,0]
# 		v_phi = m_sc_out[0,:,1]

# 		theta = m_sc_out[:,:,1].real
# 		phi = m_sc_out[:,:,0].real
# 		r = m_sc_out[:,:,6].real
			
		
# 		x = r * np.sin(phi) * np.cos(theta)
# 		y = r * np.cos(phi)
# 		z = r * np.sin(phi) * np.sin(theta)

# 		#angle to look at:
# 		phi_center = (90.0)*(np.pi/180.0)
# 		theta_center = (0.0)*(np.pi/180.0)
# 		phi_cone = 10*(np.pi/180.0)
# 		theta_cone = 10*(np.pi/180.0)
# 		indices = np.logical_and(np.absolute(phi - phi_center) < phi_cone,np.absolute(theta - theta_center) < theta_cone)
		
# 		X,Y,Z = x[indices],y[indices],z[indices]
# 		intensities.append(np.linalg.norm([X,Y,Z])**2)

# 	return intensities

if __name__ == "__main__":
	initial_particle_xyz,particle_radii = make_structure("linear_chain")
	# plot_3d_positions(particle_xyz)
	angle_min = 0
	angle_max = 360 
	angle_step = 10

	wavelengths = np.linspace(650,750,1)#[500,550,600,610] + list(np.linspace(617,617.3,5)) + [620,650,700]
	debug = 1
	scattering_at_90 = []
	scattering_at_0 = []
	angles = [0,90]#np.linspace(0,180,18)
	intensities = []
	print "Start..."
	for angle in angles:
		rotated_particle_xyz =rotate(initial_particle_xyz,0.0,0.0,angle)
		print "angle:", angle, "particles:", rotated_particle_xyz
		amplt = multiwavelength_scattering_cone(particle_xyz=rotated_particle_xyz,particle_radii=particle_radii,wavelengths=wavelengths,debug=1)
		if debug > 0:
			ax = plt.gca()
			dist = np.linalg.norm(rotated_particle_xyz[0,:]-rotated_particle_xyz[1,:])
			r1 = rotated_particle_xyz[0]
			ax.set_title("Polarization angle (rel to x): {3},Particle radius: {0} nm,\nParticle#: {1}, Separation:{2:.3g},\nPositions(x,y,z):\n".format(RADIUS,Ns,dist,gamma) + "\n".join([",".join(["{0:.3g}".format(n) for n in p]) for p in rotated_particle_xyz]))
			plt.show()	
		intensities.append(amplt)

	intensities = np.asarray(intensities)
	print "Intensity shapes:",intensities.shape
	for i in range(intensities.shape[1]):
		I = intensities[:,i]
		# I = I - np.min(I)
		# I = I/np.max(I)
		plt.plot(angles, I,'o-',color=[c/255.0 for c in wavelength_to_rgb(wavelengths[i])],label="wavelengths:{} nm".format(wavelengths[i]))
		ratio = I[1]/float(I[0])
		scattering_at_90.append(I[1])
		scattering_at_0.append(float(I[0]))
	plt.xlabel("Rotation angle (around z-axis) [deg]")
	plt.ylabel("Integrated Intensity")
	plt.title("Rotating Dimer scattering, NP radius: {0:.3g} nm, Gap Spacing: {1:.3g} nm".format(RADIUS,SPACING))
	plt.legend()
	plt.show()

	fig, axarr = plt.subplots(2,figsize=(8,16))
	axarr[0].plot(wavelengths,scattering_at_90,'o-',label="Intensity($I_{90}$),$theta=90^o$")
	axarr[0].plot(wavelengths,scattering_at_0,'o-',label="Intensity($I_{0}$),$theta=0^o$")
	ratio = np.asarray(scattering_at_90)/np.asarray(scattering_at_0)
	axarr[1].plot(wavelengths,ratio,'o-',label="$I_{90}/I_{0}$")

	axarr[0].set_xlabel("Wavelength [nm]")
	axarr[0].set_ylabel("Intensity [arb units.]")
	axarr[0].legend()
	
	axarr[1].set_xlabel("Wavelength [nm]")
	axarr[1].set_ylabel("Intensity ratio")
	axarr[1].legend()

	axarr[0].minorticks_on()
	axarr[0].grid(which="minor",linestyle="-",linewidth=0.2,color="grey")
	axarr[0].grid(which="major",linestyle="--",linewidth=0.5,color="grey")


	axarr[1].minorticks_on()
	axarr[1].grid(which="minor",linestyle="-",linewidth=0.2,color="grey")
	axarr[1].grid(which="major",linestyle="--",linewidth=0.5,color="grey")
	
	axarr[0].title("60nm AuNP dimer ratio of scattering intensities when\noriented parallel and perpendicular incident polarization direction")
	plt.show()

	# for angle in angles:
	# 	rotated_particle_xyz =rotate(initial_particle_xyz,0.0,angle,0.0)
	# 	print "scattering difference"
	# 	multiwavelength_scattering_difference(particle_xyz=rotated_particle_xyz,particle_radii=particle_radii,wl_1=wavelengths[0],wl_2=wavelengths[1])
	# plt.show()

	# intensities = [[]]*len(wavelengths)
	# def simulate_wavelength(wavelength):
	# 	return cluster_rotations(particle_xyz,angles,wavelength,particle_radii)


	# for wl in wavelengths:
	# 	simulate_wavelength(wl)

	import sys
	sys.exit()	
	# from multiprocessing import Pool 
	# pool = Pool()
	# intensities = pool.map(simulate_wavelength,wavelengths)

	for i in range(0,len(wavelengths)):
		wl = wavelengths[i]
		color = wavelength_to_rgb(wl)
		intensity = np.asarray(intensities[i])
		intensity = intensity - np.min(intensity)
		intensity = intensity/np.max(intensity)
		plt.plot(angles,intensity,'-x',color=[c/255.0 for c in color],label="Wavelength:{} nm".format(wl))
	plt.title("AuNP Random, NP radius: {0} nm, NPs: {1}".format(RADIUS,Ns))
	plt.xlabel("Rotation angle [deg]")
	plt.ylabel("Normalized intensity (range [0:1]) [arb. units]")
	plt.legend()
	plt.show()
