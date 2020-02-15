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
sigma1 = np.array([[1, 0], [0, 1]])
mu = np.array([0, 0])
mu1 = np.array([2, 2])
mu2 = np.array([2, -2])
mu3 = np.array([-2, 2])
mu4 = np.array([-2, -2])

# Variable variance.
sigma = np.array([[10, 0], [0, 10]])

# List means.
means = [mu1, mu2, mu3, mu4]

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
Z = multivariate_gaussian(pos, mu, sigma)

for mean in means:
    Z += multivariate_gaussian(pos, mean, sigma1)

# The distribution on the variables X, Y packed into pos.
ax.plot_surface(X, Y, Z, rstride=3, cstride=3, linewidth=1,
                antialiased=False, cmap=cm.viridis, zorder=0)
cset = ax.contourf(X, Y, Z, zdir='z', offset=-0.3, cmap=cm.viridis)

# Adjust the limits, ticks and view angle
ax.set_zlim(-0.3, 0.4)
ax.set_zticks(np.linspace(0, 0.4, 5))
ax.view_init(27, -21)

plt.show()
