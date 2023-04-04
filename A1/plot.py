# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 
import json

# %%
print("Reading in data...")

# read in the constants
with open("data/output_constants.json", "r") as f:
    consts = json.load(f)
M = consts["M"]
x_min = consts["x_min"]
x_max = consts["x_max"]
h = consts["dx"]

# read in the simulation data
data = np.genfromtxt("data/output.csv", delimiter=",")
n = data[:,0]
t = data[:,1]
E = data[:,2:2+M]
B = data[:,2+M:2+2*M]

x_E = np.linspace(x_min,       x_min + (M-1.0)*h, M)
x_B = np.linspace(x_min+0.5*h, x_min + (M-0.5)*h, M)

# %%
print("Plotting...")
fig = plt.figure() 
ax = plt.axes(xlim=(x_min, x_max), ylim=(-1.5, 1.5))
ax.set_xlabel("x")
line_E, = ax.plot([], [], lw=2, label="E")
line_B, = ax.plot([], [], lw=2, label="B")
ax.legend(loc="upper right")

frame_count = E.shape[0]

def init():  
    line_E.set_data([], [])
    line_B.set_data([], []) 
    return line_E, line_B 

def animate(frame_idx):
    if (frame_idx % 10 == 0) or (frame_idx == frame_count-1): 
        print(f"\rRendering... ({100*((frame_idx+1)/frame_count):.1f}%)", end="")
    ax.set_title(f"t = {t[frame_idx]:.2f}")
    line_E.set_data(x_E, E[frame_idx,:])
    line_B.set_data(x_E, B[frame_idx,:])
    return line_E, line_B

anim = animation.FuncAnimation(
	fig, animate, init_func=init, 
	frames=frame_count, interval=20, blit=True) 

anim.save('figures/movie.gif', writer='imagemagick')
print("\nRendering done!")
