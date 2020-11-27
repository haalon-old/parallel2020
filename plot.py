from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import pandas as pd


df = pd.read_csv('data-64-25.csv')
# df.set_index(['x', 'y', 'z'])
df2 = df.loc[df['z'] == 24].to_numpy()

df2 = df2[df2[:,1].argsort(kind='mergesort')]
df2 = df2[df2[:,0].argsort(kind='mergesort')]
print(df2[:9])

xs = (df2[:,0] / 64).reshape((64,64))
ys = (df2[:,1] / 64).reshape((64,64))
zs = df2[:,3].reshape((64,64))
# xs, ys = np.meshgrid(xs, ys)
print(xs.shape)
print(len(xs))



t = 0.25
n = 64

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(0, 64, 1)
Y = np.arange(0, 64, 1)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)

# Plot the surface.
surf = ax.plot_surface(xs, ys, zs, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

