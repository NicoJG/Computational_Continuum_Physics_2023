# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 
import json

# %%
print("Reading in data...")
E_arr = np.genfromtxt("data/output_E.csv")
B_arr = np.genfromtxt("data/output_B.csv")

# read in the constants
with open("data/output_constants.json", "r") as f:
    consts = json.load(f)

M = consts["M"]
x_min = consts["x_min"]
x_max = consts["x_max"]
h = consts["dx"]

N = consts["N"]
t_min = consts["t_min"]
t_max = consts["t_max"]
tau = consts["dt"]

assert N == E_arr.shape[0] == B_arr.shape[0]
assert M == E_arr.shape[1] == B_arr.shape[1]

x_E = np.linspace(x_min,       x_min + (M-1.0)*h, M)
x_B = np.linspace(x_min+0.5*h, x_min + (M-0.5)*h, M)

# %%
print("Plotting...")
fig = plt.figure() 
ax = plt.axes(xlim=(x_min, x_max), ylim=(-2, 2))
ax.set_xlabel("x")
line_E, = ax.plot([], [], lw=2, label="E")
line_B, = ax.plot([], [], lw=2, label="B")
ax.legend(loc="upper right")

frame_count = 200
frame_stride  = N//frame_count
assert(frame_count * frame_stride <= N)

def init():  
    line_E.set_data([], [])
    line_B.set_data([], []) 
    return line_E, line_B 

def animate(n):
    if (n % 10 == 0): print(f"\rRendering... ({100*((n+10)/frame_count):.1f}%)", end="")
    ax.set_title(f"t = {tau*n*frame_stride:.2f}")
    line_E.set_data(x_E, E_arr[n*frame_stride,:])
    line_B.set_data(x_E, B_arr[n*frame_stride,:])
    return line_E, line_B

anim = animation.FuncAnimation(
	fig, animate, init_func=init, 
	frames=frame_count, interval=20, blit=True) 

anim.save('figures/movie.gif', writer='imagemagick')
print("\nRendering done!")
