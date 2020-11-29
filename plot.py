from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import pandas as pd
import math

rcParams['savefig.pad_inches'] = 0


t = 0.25
n = 64
l_x = 1
l_y = 1
l_z = 1

def u(x, y, z):
    a_t = math.pi * math.sqrt(1/(l_x*l_x) + 4/(l_y*l_y) + 9/(l_z*l_z))
    return np.sin(math.pi * x / l_x) * np.sin(2 * math.pi * y / l_x) * np.sin(3 * math.pi * z / l_x) * math.cos(t * a_t)


def plot(ax1, ax2, ax3, sx, sy, i, title=""):
    ax = fig.add_subplot(sx ,sy ,i, projection='3d')
    ax.set_title(title)
    # Plot the surface.
    surf = ax.plot_surface(ax1, ax2, ax3, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    # # # Customize the z axis.
    # # ax.set_zlim(-1, 1)
    # # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))



df = np.genfromtxt('data-64-25.csv', delimiter=',')
# df.set_index(['x', 'y', 'z'])


df = df[df[:,0].argsort(kind='mergesort')]
df = df[df[:,1].argsort(kind='mergesort')]
df = df[df[:,2].argsort(kind='mergesort')]

df[:,0] /= n
df[:,1] /= n
df[:,2] /= n

print(df.shape)

true_z = u(df[:,0], df[:,1], df[:,2])

df = np.hstack((df, np.atleast_2d(true_z).T))

err = np.abs(df[:,3] - df[:,4])
df = np.hstack((df, np.atleast_2d(err).T))

print(df[-4:])

max_i = np.argmax(err)
print(max_i)
print(df[max_i])

max_x, max_y, max_z = df[max_i, :3]
print(max_x, max_y, max_z)

fig = plt.figure()

df2 = df[np.where(df[:,0] == max_x)]

ys = df2[:,1].reshape((64,64))
zs = df2[:,2].reshape((64,64))
tvs = df2[:,4].reshape((64,64))
avs = df2[:,3].reshape((64,64))
es = df2[:,5].reshape((64,64))

plot(ys,zs,avs,3,3,1, f"Approximated plot at x={max_x}")
plot(ys,zs,tvs,3,3,2, f"True plot at x={max_x}")
plot(ys,zs,es,3,3,3, f"Error plot at x={max_x}")


df2 = df[np.where(df[:,1] == max_y)]

xs = df2[:,0].reshape((64,64))
zs = df2[:,2].reshape((64,64))
tvs = df2[:,4].reshape((64,64))
avs = df2[:,3].reshape((64,64))
es = df2[:,5].reshape((64,64))

plot(xs,zs,avs,3,3,4, f"Approximated plot at y={max_y}")
plot(xs,zs,tvs,3,3,5, f"True plot at y={max_y}")
plot(xs,zs,es,3,3,6, f"Error plot at y={max_y}")

df2 = df[np.where(df[:,2] == max_z)]

ys = df2[:,1].reshape((64,64))
zs = df2[:,2].reshape((64,64))
tvs = df2[:,4].reshape((64,64))
avs = df2[:,3].reshape((64,64))
es = df2[:,5].reshape((64,64))

plot(xs,ys,avs,3,3,7, f"Approximated plot at z={max_z}")
plot(xs,ys,tvs,3,3,8, f"True plot at z={max_z}")
plot(xs,ys,es,3,3,9, f"Error plot at z={max_z}")



# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)

plt.autoscale(tight=True)

plt.show()

