import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation 

# Potential to plot level curves for
def V(x, y):
    return -5 / (1 + (x/5)**2 + (y/4)**2)**4

# Folder with csv files
listing = sorted(os.listdir("data"))
listing.remove("psi_time_and_frequency.csv")
psi0 = np.genfromtxt("data/" + listing[0], delimiter=",")
M = psi0.shape[0]
L = 10

x = ((np.array(range(0,M))-0.5)*(2.0/M)-1.0)*L
y = ((np.array(range(0,M))-0.5)*(2.0/M)-1.0)*L
X, Y = np.meshgrid(x, y)

clim = 0.5
cmap = "seismic"

# Save each frame separately
for (i, file) in enumerate(listing):
    print(f"Generating figure {i+1}/{len(listing)}")
    psi = np.genfromtxt("data/" + file, delimiter=",")

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect("equal")
    im = ax.pcolormesh(X, Y, psi, vmin=-clim, vmax=clim, cmap=cmap, rasterized=True)
    ax.contour(X, Y, V(X,Y), colors="black", linestyles="solid", alpha=0.2, levels=[-4, -3, -2, -1])
    fig.colorbar(im, ax=ax)
    fig.savefig("figures/" + file.split(".")[0] + "." + file.split(".")[1]  + ".pdf")

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
    print(f"Rendering frame {frame_idx+1}/{len(listing)}")
    file = listing[frame_idx]
    psi = np.genfromtxt("data/" + file, delimiter=",")
    im = ax.pcolormesh(X, Y, psi, vmin=-clim, vmax=clim, cmap=cmap)
    return im

anim = animation.FuncAnimation(
	fig, animate, init_func=init, 
	frames=len(listing), interval=100)

anim.save('figures/movie.gif', writer='imagemagick')
