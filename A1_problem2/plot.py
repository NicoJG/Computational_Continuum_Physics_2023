import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.animation as animation 

# Folder with csv files
listing = sorted(os.listdir("data"))
B0 = np.genfromtxt("data/" + listing[0], delimiter=",")
M = B0.shape[0]

clim = 1
cmap = "seismic"

# Create animation
fig = plt.figure() 
ax = fig.add_subplot(1,1,1)
ax.set_aspect("equal")
im = ax.pcolormesh(B0, vmin=-clim, vmax=clim, cmap=cmap)
fig.colorbar(im, ax=ax)

def init():  
    im = ax.pcolormesh(B0, vmin=-clim, vmax=clim, cmap=cmap)
    return im

def animate(frame_idx):
    file = listing[frame_idx+1]
    B = np.genfromtxt("data/" + file, delimiter=",")
    im = ax.pcolormesh(B, vmin=-clim, vmax=clim, cmap=cmap)
    return im

anim = animation.FuncAnimation(
	fig, animate, init_func=init, 
	frames=len(listing)-1, interval=500)

anim.save('figures/movie.gif', writer='imagemagick')


# Save each frame separately
for file in listing:
    print(file)
    B = np.genfromtxt("data/" + file, delimiter=",")

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect("equal")
    im = ax.pcolormesh(B, vmin=-clim, vmax=clim, cmap=cmap)
    fig.colorbar(im, ax=ax)
    fig.savefig("figures/" + file.split(".")[0] + "." + file.split(".")[1]  + ".png")

    plt.close(fig)

