import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import json
from pathlib import Path
import shutil
from tqdm import tqdm

animation_dir = Path("figures/1d_animation")
if animation_dir.is_dir():
    shutil.rmtree(animation_dir)
animation_dir.mkdir(parents=True)

# read in the constants
with open("data/1d_simulation.csv", "r") as f:
    consts = json.loads(f.readline())

x_min = consts["x_min"]
x_max = consts["x_max"]
M = consts["M"]
x = np.linspace(x_min, x_max, M)

data = np.genfromtxt("data/1d_simulation.csv", delimiter=",", skip_header=2)

t = data[:, 0]
u = data[:,1:]

for frame in [0,len(t)-1]:
    plt.figure()
    plt.plot(x,u[frame,:])
    plt.tight_layout()
    plt.savefig(animation_dir/f"frame_{frame:04d}.png")
    plt.close()

# Create animation
frames = len(t)
fps = 10

fig,ax = plt.subplots()
progress = tqdm(total=frames)
def animate(frame):
    progress.update()
    ax.clear()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(u[frame].min(),u[frame].max())
    ax.set_title(f"t = {t[frame]}")
    line, = ax.plot(x,u[frame,:])
    return line,
        
ani = FuncAnimation(fig, animate, interval=1000/fps, blit=True, repeat=True, frames=frames)    
ani.save("figures/1d_simulation.gif", dpi=300, writer=PillowWriter(fps=fps))

