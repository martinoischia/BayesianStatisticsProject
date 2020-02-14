import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


# Code base from https://scipython.com/blog/
# visualizing-the-bivariate-gaussian-distribution/.

N = 60
X = np.linspace(-10, 10, N)
Y = np.linspace(-10, 10, N)
X, Y = np.meshgrid(X, Y)

# Fixed Distributions parameters.
sigma = np.array([[1, 0], [0, 1]])
mu1 = np.array([0, 0])
mu2 = np.array([-3, 3])
mu3 = np.array([3, -3])

# Variable mean.
mu = np.array([5, 5])

# List means.
means = [mu, mu1, mu2, mu3]

# Pack X and Y into a single 3-dimensional array
pos = np.empty(X.shape + (2,))
pos[:, :, 0] = X
pos[:, :, 1] = Y


def multivariate_gaussian(pos, mu, sigma):

    n = mu.shape[0]
    sigma_det = np.linalg.det(sigma)
    sigma_inv = np.linalg.inv(sigma)
    num = np.sqrt((2*np.pi)**n * sigma_det)
    # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
    # way across all the input variables.
    fac = np.einsum('...k,kl,...l->...', pos-mu, sigma_inv, pos-mu)

    return np.exp(-fac / 2) / num


fig = plt.figure()
ax = fig.gca(projection='3d')
Z = 0

for mean in means:
    Z += multivariate_gaussian(pos, mean, sigma)

# The distribution on the variables X, Y packed into pos.
ax.plot_surface(X, Y, Z, rstride=3, cstride=3, linewidth=1,
                antialiased=False, cmap=cm.viridis, zorder=0)
cset = ax.contourf(X, Y, Z, zdir='z', offset=-0.15, cmap=cm.viridis)

# Adjust the limits, ticks and view angle
ax.set_zlim(-0.15, 0.2)
ax.set_zticks(np.linspace(0, 0.2, 5))
ax.view_init(27, -21)

plt.show()
