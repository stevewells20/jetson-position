"""
============
3D animation
============

Units of Measure:
	Force => (Kilogram * Meter) / Second^2
	Length => Meter
	Mass => Kilogram
	Time => Second
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import cm

import time

class AxialForce(object):
	"""docstring for AxialForce."""
	def __init__(self, arg):
		super(AxialForce, self).__init__()
		self.arg = arg


class ROV(object):
	"""
	docstring for ROV.
	"""

	def __init__(self, lims=2.5):
		super(ROV, self).__init__()
		# self.start_time = time.perf_counter()
		self.bounds = [lims, lims, lims]
		self.init_vars()
		# self.stop_time = time.perf_counter()

	def init_vars(self):
		'''
		p => meters
		v => meters / second
		a => meters / second^2
		'''
		self.p = np.array([0, 0, 0])  # [x,y,z]
		self.v = np.array([0, 0, 10])  # [x,y,z]
		self.a = np.array([0, 0, 0])  # [x,y,z]
		self.mass = 1.0  # kg
		self.forces = np.array([0, 0, 0])  # [x,y,z]
		self.force_line, = ax.plot([0,0], [0,0], zs=[0,0])
		self.vel_line, = ax.plot([0,0], [0,0], zs=[0,0])
		self.drag_line, = ax.plot([0,0], [0,0], zs=[0,0])
		self.target_p, = ax.plot([self.p[0]],[self.p[1]],[self.p[2]], 'h')
		self.current_p, = ax.plot([self.p[0]], [self.p[1]], [self.p[2]], 'H')
		self.x_line, = ax.plot([-self.bounds[0],self.p[0]] , \
							   [self.p[1],self.p[1]], \
							   zs=[self.p[2],self.p[2]])
		self.y_line, = ax.plot([self.p[0],self.p[0]], \
							   [self.bounds[1],self.p[1]], \
							   zs=[self.p[2],self.p[2]])
		self.z_line, = ax.plot([self.p[0],self.p[0]], \
							   [self.p[1],self.p[1]], \
							   zs=[-self.bounds[2],self.p[2]])
		self.dt = 0.1
		self.elasticity_coeff = 1.0

	def drag_force(self, v_magnitude):
		'''
		F_D => drag_force
		Cd => drag_coeff (sphere=0.47,
							half_sphere=0.42,
							l_cylinder=0.82,
							model_rocket=0.75)
		rho => density (water= 1000kg / m^3) [change rho for diff environs]
		V => flow_velocity (self.v possibly?): our vel in relation to water
		A => reference area
		F_D = Cd * A * 1/2( rho * V^2 )

		'''
		drag_coeff = 0.82 # dimensionless coeff, 0.82 for l_cylinder
		rho = 1000.0 # kg / m^3
		# if looking at ROV from front, is visible surface area
		ref_area = 	1.0 # m^2 , projected frontal area of the vehicle
		force_drag = drag_coeff * ref_area * 0.5 * ( rho * (v_magnitude ** 2) )
		# print ('force_drag: ', force_drag)
		return force_drag

	def reset_vars(self):
		self.p = np.array([0, 0, 0])  # [x,y,z]
		self.v = np.array([0, 0, 0])  # [x,y,z]
		self.a = np.array([0, 0, 0])  # [x,y,z]
		self.forces = np.array([0, 0, 0])  # [x,y,z]

	def time_start(self):
		self.start_time = time.perf_counter()

	def time_stop(self):
		self.stop_time = time.perf_counter()
		return (self.stop_time - self.start_time)

	def sph2cart(self, sph_coord):
		az, el, r = sph_coord
		rcos_theta = r * np.cos(el)
		x = rcos_theta * np.cos(az)
		y = rcos_theta * np.sin(az)
		z = r * np.sin(el)
		# print (sph_coord,(x,y,z))
		return [x, y, z]

	def mag_cap(self, val, max_magnitude=2.5):
		sign = 1
		if val < 0:
			sign = -1
		out = sign* max(min(abs(val), max_magnitude),0)
		# print ("(val, out) : ",val, out)
		return out

	def cart2sph(self, cart_coord):
		x, y, z = cart_coord
		hxy = np.hypot(x, y)
		r = np.hypot(hxy, z)
		el = np.arctan2(z, hxy)
		az = np.arctan2(y, x)
		sph_coord = [az,el,r]
		# print("cart2sph: ",cart_coord, sph_coord)
		return [az, el, r]

	def update_loc(self, n):
		v_magnitude = np.linalg.norm(self.v) # m / s
		print ("v_magnitude: ",v_magnitude)

		v_sph = self.cart2sph(self.v) # x,y,z => az,el,r
		drag_sph = v_sph
		drag_sph[2] = self.drag_force(v_magnitude)
		drag_cart = self.sph2cart(drag_sph)
		self.drag_line.set_data([self.p[0],drag_cart[0]],[self.p[1],drag_cart[1]])
		self.drag_line.set_3d_properties([self.p[2],drag_cart[2]], 'z')

		self.a = ((self.forces / self.mass)*self.dt)
		self.v = self.v + self.a*self.dt

		for axis in range(len(self.p)):
			if abs(self.p[axis]) >= self.bounds[axis]:
				self.v[axis] = -self.v[axis]*self.elasticity_coeff

		self.p = self.p + self.v*self.dt

		# print ("x=",self.a[0], "\ty=", self.a[1], "\tz=", self.a[2], "\t")
		self.current_p.set_data(np.array([self.p[0], self.p[1]]))
		self.current_p.set_3d_properties(self.p[2], 'z')

		self.x_line.set_data([-self.bounds[0],self.p[0]],[self.p[1],self.p[1]])
		self.x_line.set_3d_properties([self.p[2],self.p[2]], 'z')
		self.y_line.set_data([self.p[0],self.p[0]], [self.bounds[1],self.p[1]])
		self.y_line.set_3d_properties([self.p[2],self.p[2]], 'z')
		self.z_line.set_data([self.p[0],self.p[0]],[self.p[1],self.p[1]])
		self.z_line.set_3d_properties([-self.bounds[2],self.p[2]], 'z')

		return [self.current_p, self.x_line, self.y_line]

	def press(self, event):
		print('press', event.key)
		# sys.stdout.flush()
		if event.key == '7':
			self.forces[0] = 10
		if event.key == '4':
			self.forces[0] = 0
		if event.key == '1':
			self.forces[0] = -10
		if event.key == '8':
			self.forces[1] = 10
		if event.key == '5':
			self.forces[1] = 0
		if event.key == '2':
			self.forces[1] = -10
		if event.key == '9':
			self.forces[2] = 10
		if event.key == '6':
			self.forces[2] = 0
		if event.key == '3':
			self.forces[2] = -10
		if event.key == '/':
			self.reset_vars()

# def update_lines(n, point):
#     point.set_data(np.array([x[n], y[n]]))
#     point.set_3d_properties(z[n], 'z')
#     return point

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)
ax_dev = 2.5
rov = ROV(ax_dev);

# Connect Callbacks
fig.canvas.mpl_connect('key_press_event', rov.press)

# Setting the axes properties

ax.set_xlim3d([-ax_dev, ax_dev])
ax.set_ylim3d([-ax_dev, ax_dev])
ax.set_zlim3d([-ax_dev, ax_dev])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Test')

# Creating the Animation object
# point, = ax.plot([rov.p[0]], [rov.p[1]], [rov.p[2]], 'o')
ani = animation.FuncAnimation(fig, rov.update_loc, 99,			\
								interval=10, blit=False)

plt.show()
