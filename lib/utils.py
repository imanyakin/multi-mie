import numpy as np 

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