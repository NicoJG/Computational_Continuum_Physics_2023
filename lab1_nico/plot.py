# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 

# %%
# read data
u = np.genfromtxt("data/output.csv", delimiter=",")

# constants of the simulation
N = u.shape[0]
M = u.shape[1]
x_min = 0
x_max = 1
dx = (x_max - x_min) / (M-1)
t_min = 0
t_max = 10
dt = (t_max - t_min) / (N-1)

x = np.linspace(x_min, x_max, M)

# %%
print("Plotting...")
fig = plt.figure() 
ax = plt.axes(xlim=(x_min, x_max), ylim=(0, 1))
ax.set_xlabel("x")
ax.set_ylabel("u")
line, = ax.plot([], [], lw=2)

frame_count = 200
frame_stride  = N//frame_count
assert(frame_count * frame_stride <= N)

def init():  
    line.set_data([], [])
    return (line,)

def animate(n):
    if (n % 10 == 0): print(f"\rRendering... ({100*((n+10)/frame_count):.1f}%)", end="")
    ax.set_title(f"t = {dt*n*frame_stride:.2f}")
    line.set_data(x, u[n*frame_stride,:])
    return (line,)

anim = animation.FuncAnimation(
	fig, animate, init_func=init, 
	frames=frame_count, interval=20, blit=True) 

anim.save('figures/movie.gif', writer='imagemagick')
print("\nRendering done!")
