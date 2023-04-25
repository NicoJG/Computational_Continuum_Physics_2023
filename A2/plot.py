import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation 


def V(x, y):
    return -5 / (1 + (x/5)**2 + (y/4)**2)**4

# Folder with csv files
listing = sorted(os.listdir("data"))
psi0 = np.genfromtxt("data/" + listing[0], delimiter=",")
M = psi0.shape[0]

x = np.linspace(-10, 10, M)
y = np.linspace(-10, 10, M)
X, Y = np.meshgrid(x, y)

clim = 1
cmap = "seismic"

# Save each frame separately
for file in listing:
    print(file)
    psi = np.genfromtxt("data/" + file, delimiter=",")

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect("equal")
    im = ax.pcolormesh(X, Y, psi, vmin=-clim, vmax=clim, cmap=cmap)
    ax.contour(X, Y, V(X,Y), colors="black", linestyles="solid", alpha=0.2, levels=[-4, -3, -2, -1])
    fig.colorbar(im, ax=ax)
    fig.savefig("figures/" + file.split(".")[0] + "." + file.split(".")[1]  + ".png")

    plt.close(fig)


# Create animation
fig = plt.figure() 
ax = fig.add_subplot(1,1,1)
ax.set_aspect("equal")
ax.contour(X, Y, V(X,Y), colors="black", linestyles="solid", alpha=0.2, levels=[-4, -3, -2, -1])
fig.colorbar(im, ax=ax)

def init():  
    im = ax.pcolormesh(X, Y, psi0, vmin=-clim, vmax=clim, cmap=cmap)
    return im

def animate(frame_idx):
    file = listing[frame_idx]
    psi = np.genfromtxt("data/" + file, delimiter=",")
    im = ax.pcolormesh(X, Y, psi, vmin=-clim, vmax=clim, cmap=cmap)
    return im

anim = animation.FuncAnimation(
	fig, animate, init_func=init, 
	frames=len(listing), interval=500)

anim.save('figures/movie.gif', writer='imagemagick')


