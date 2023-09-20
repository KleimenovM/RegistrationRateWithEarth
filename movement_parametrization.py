# Neutrino through Earth propagation
# Parametrize the source daily movement

import numpy as np
import matplotlib.pyplot as plt

from telescope import get_simple_telescope
from source import get_sources
from tools import deg_to_rad, sph_coord, rot_matrix


def plot_a_sphere(axis: plt.axis, rotation_angle: float):
	"""
	Creates a sphere with parallels and meridians on it
	@param axis: an axis on which the sphere is to be drawn
	@param rotation_angle: an angle (to vertical direction) at which the sphere is to be turned
	@return:
	"""
	# meridians
	shift = deg_to_rad([20])
	cur_angle = 0
	n = 100
	theta = np.linspace(-np.pi/2, np.pi/2 + 2 * np.pi / n, n)

	rm = rot_matrix(rotation_angle)

	while cur_angle < 2 * np.pi:
		vec = np.zeros([3, n])
		vec[0], vec[1], vec[2] = sph_coord(r=1, theta=theta, phi=cur_angle)

		vec = np.dot(rm, vec)

		axis.plot(vec[0], vec[1], vec[2], color='lightgray')
		cur_angle += shift

	# parallels
	shift = deg_to_rad([15])
	cur_angle = -np.pi / 2
	phi = np.linspace(0, 2 * np.pi, n)

	while cur_angle < np.pi / 2:
		vec = np.zeros([3, n])

		vec[0], vec[1], vec[2] = sph_coord(r=1, theta=cur_angle, phi=phi)

		vec = np.dot(rm, vec)

		if abs(cur_angle) < 1e-3:
			axis.plot(vec[0], vec[1], vec[2], color='black')
		else:
			axis.plot(vec[0], vec[1], vec[2], color='lightgray')

		cur_angle += shift

	m1 = 10
	vec_vert = np.zeros([3, m1])
	vec_vert[2] = np.linspace(-1.1, 1.1, m1)
	axis.plot(vec_vert[0], vec_vert[1], vec_vert[2], color='red', linewidth=2, linestyle='dashed')

	vec = np.zeros([3, n])
	vec[0], vec[1], vec[2] = sph_coord(r=1, theta=0, phi=np.linspace(0, 2*np.pi, 100))
	axis.plot(vec[0], vec[1], vec[2], color='red')

	vec_vert2 = np.dot(rm, vec_vert)
	axis.plot(vec_vert2[0], vec_vert2[1], vec_vert2[2], color='black', linewidth=2, linestyle='dashed')

	return
		

if __name__ == '__main__':
	baikal_latitude = deg_to_rad([51, 46])  # 51 46' N to rad
	source_numbers = [0, 1, 2, 3, 5, 6, 9, 10]

	sources = get_sources("data/source_table.csv")

	t = get_simple_telescope("Baikal", [51, 46], "data/eff_area.root")

	fig = plt.figure(figsize=(12, 6))

	ax = fig.add_subplot(121, projection='3d')
	ax2 = fig.add_subplot(122)
	plot_a_sphere(ax, baikal_latitude)

	ax.axis('scaled')
	ax.axis('off')

	ax2.plot((0, 2 * np.pi), (0, 0), color='black')
	ax2.set_xlabel(r'$\psi,\ rad$')
	ax2.set_ylabel(r'$\theta,\ rad$')

	psi_sample = np.linspace(0, 2 * np.pi, 180)

	for n in source_numbers:
		s_i = sources[n]
		print(s_i.info())

		vec1, theta1 = t.get_orbit_parametrization(source=s_i, angular_precision=180)

		ax.plot(vec1[0], vec1[1], vec1[2], label=s_i.name)

		ax2.plot(psi_sample, theta1)

	ax2.plot(psi_sample, np.ones(psi_sample.size) * (-0.5))
	ax.legend()
	plt.tight_layout()
	# plt.savefig("vis.png")
	plt.show()
