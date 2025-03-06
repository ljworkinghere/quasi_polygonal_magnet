import numpy as np
from scipy import integrate
from scipy.integrate import quad
from scipy import optimize
from scipy.interpolate import interp1d
import time
from matplotlib.path import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# ax = plt.figure().add_subplot(111, projection='3d')


# generate the CCT path that is used in project-rat 


# calculate the area of the CCT former 
def polygon_area(vertices):
    n = len(vertices)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        xi, yi = vertices[i]
        xj, yj = vertices[j]
        area += (xi * yj) - (xj * yi)
    return abs(area) / 2.0



# inner CCT former
rho_1 = 0.6
N = int(1e4)
theta_array = 2*np.pi*np.arange(N)/N

vertices_1 = []
for i in range(N):
    z = rho_1 * np.exp(1j * theta_array[i])
    z_eta = 0.5*(z + 1/z) # conformal mapping
    vertices_1.append([z_eta.real, z_eta.imag])
vertices_1 = np.array(vertices_1)

area_total = polygon_area(vertices_1)
print("Total area is:", area_total)
plt.plot(vertices_1[:, 0], vertices_1[:, 1], color="r")
plt.axis("equal")
plt.show()

pitch = 0.12 # the pitch of the CCT dipole coil
Amp1 = 0.9 # the amplitude of the CCT dipole coil
theta_array = np.linspace(-100*np.pi, 100*np.pi, 1+40*2000) # 100 turns
res = []
scalar_coeff_1 = 1/np.sqrt(area_total) # normalize the area
for theta in theta_array:
    z = rho_1 * np.exp(1j * theta)
    z_eta = 0.5*(z + 1/z)*scalar_coeff_1 # conformal mapping
    x_cor = z_eta.real
    y_cor = z_eta.imag
    z_cor = np.sin(1*theta)*Amp1 + pitch*(theta/2/np.pi)
    xyz_cor = np.array([x_cor, y_cor, z_cor])
    direction = np.array([-0.5*rho_1*np.sin(1*theta)*scalar_coeff_1-0.5*np.sin(theta)/rho_1*scalar_coeff_1,\
                          0.5*rho_1*np.cos(1*theta)*scalar_coeff_1-0.5*np.cos(theta)/rho_1*scalar_coeff_1,\
                          1*np.cos(1*theta)*Amp1+pitch*(1/2/np.pi)])
    direction = direction/np.linalg.norm(direction)
    norm_vec1 = [-0.5*rho_1*np.cos(1*theta)+0.5*np.cos(theta)/rho_1,\
                -0.5*rho_1*np.sin(1*theta)-0.5*np.sin(theta)/rho_1,\
                0]
    norm_vec1 = norm_vec1 / np.linalg.norm(norm_vec1)
    res.append(np.hstack((xyz_cor, direction, norm_vec1)))
    # In project-rat, the path needs 3*3 Coordinates: 
    # Central Coordinates of Current Element;
    # Its Direction Vector;
    # One Normal Vector.
res = np.array(res)
# from mpl_toolkits.mplot3d import Axes3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(res[:, 0],res[:, 1], res[:, 2], c='r', marker='o')
# plt.show()
file_name = "./Elliptical_cross_one_turn/" +"level1.txt"
np.savetxt(file_name, res, delimiter="    ")



# outer CCT former, repeate the above precedure
pitch = 0.12
Amp1 = -0.9
theta_array = np.linspace(-100*np.pi, 100*np.pi, 1+40*2000)
res = []
scalar_coeff_1 = 1.05/np.sqrt(area_total) # normalize the former area to 1.05^2
for theta in theta_array:
    z = rho_1 * np.exp(1j * theta)
    z_eta = 0.5*(z + 1/z)*scalar_coeff_1 # conformal mapping
    x_cor = z_eta.real
    y_cor = z_eta.imag
    z_cor = np.sin(1*theta)*Amp1 + pitch*(theta/2/np.pi)
    xyz_cor = np.array([x_cor, y_cor, z_cor])
    direction = np.array([-0.5*rho_1*np.sin(1*theta)*scalar_coeff_1-0.5*np.sin(theta)/rho_1*scalar_coeff_1,\
                          0.5*rho_1*np.cos(1*theta)*scalar_coeff_1-0.5*np.cos(theta)/rho_1*scalar_coeff_1,\
                          1*np.cos(1*theta)*Amp1+pitch*(1/2/np.pi)])
    direction = direction/np.linalg.norm(direction)
    norm_vec1 = [-0.5*rho_1*np.cos(1*theta)+0.5*np.cos(theta)/rho_1,\
                -0.5*rho_1*np.sin(1*theta)-0.5*np.sin(theta)/rho_1,\
                0]
    norm_vec1 = norm_vec1 / np.linalg.norm(norm_vec1)
    res.append(np.hstack((xyz_cor, direction, norm_vec1)))
    # In project-rat, the path needs 3*3 Coordinates: 
    # Central Coordinates of Current Element;
    # Its Direction Vector;
    # One Normal Vector.
res = np.array(res)
file_name = "./Elliptical_cross_one_turn/" +"level2.txt"
np.savetxt(file_name, res, delimiter="    ")
