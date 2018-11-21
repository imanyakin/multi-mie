import numpy as np 

RADIUS  = 30
SPACING = 1
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

def generate_linear_chain(start,stop,radius,spacing):
	#generate vector in direction of chain
	direction = stop-start
	direction = direction/np.linalg.norm(direction)

	#length
	length = np.linalg.norm(stop-start)

	if length < 2*radius:
		raise ValueError("Segment length less than one particle diameter - cant fit particle into segment")
	else:
		particle_number = int(np.floor(length/float(2*radius + spacing)))
		positions = []
		for i in range(1,particle_number+1):
			if i == 1:
				pos = start + direction*radius
			else:
				pos = positions[-1]+ direction*2*radius + direction*spacing
			
			positions.append(pos)

		return np.array(positions)




def structure_generator(case="dimer", radius=RADIUS,spacing=SPACING,Ns = 2): #radius and spacing in nm
	
	if case == "dimer":
		start = np.array([0,0,0])
		stop = 4*radius + 1*spacing * np.array([1,0,0])
		xyz = generate_linear_chain(start,stop,radius,spacing)
		radii = np.array([radius]*2)
		xyz = center_COM(xyz)

	elif case == "knight":
		chain0_start = np.array([0,0,0])
		chain0_stop = chain0_start + (8*radius + 3*spacing)*np.array([1,0,0])
		chain0 = generate_linear_chain(chain0_start,chain0_stop,radius,spacing)
		
		chain1_start = chain0[-1] + np.array([0,radius + spacing,0])
		chain1_stop = chain1_start + (3*radius + spacing)*np.array([0,1,0]) 
		
		chain1 = generate_linear_chain(chain1_start,chain1_stop, radius,spacing)
		
		xyz = np.concatenate([chain0,chain1])
		xyz = center_COM(xyz)
		radii = np.array([radius]*(len(chain0)+len(chain1)))
	
	return xyz,radii

def rotate(points, theta_x,theta_y, theta_z):
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

	outp = rotate_x(points)
	outp = rotate_y(outp)
	outp = rotate_z(outp)
	return outp

if __name__ == "__main__":
	from plot_functions import plot_3d_positions
	import matplotlib.pyplot as plt 
	radius = 30 
	spacing = 1
	knight_xyz,knight_rs = structure_generator("knight")
	print knight_xyz.shape
	print knight_rs.shape
	angles = range(0,370,10)
	def plot_angle(angle):
		fig = plt.figure(figsize=(8,8))
		ax = fig.add_subplot(111, projection='3d')
		xyz = rotate(knight_xyz,angle,0,0)
		plot_3d_positions(xyz,knight_rs,ax=ax,show_plot=False,with_mayavi=False)
		ax.set_xlim(-150,150)
		ax.set_ylim(-150,150)
		ax.set_zlim(-150,150)
		plt.savefig("../images/tmp/X_knight{}.png".format(angle))
		plt.close()
        
	for angle in angles:
		print angle
		plot_angle(angle)

