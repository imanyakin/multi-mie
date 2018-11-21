import numpy as np 

from structures import structure_generator, rotate

initial_particle_xyz,particle_radii = structure_generator("dimer")
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
