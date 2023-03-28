# %%
import numpy as np
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

# %%
E_arr = np.genfromtxt("output_E.csv")
B_arr = np.genfromtxt("output_B.csv")

x_min = -1
x_max = 1
N = E_arr.shape[0]
M = E_arr.shape[1]

h = (x_max-x_min)/(M-0.5)
tau = 0.5*h

x_E = np.linspace(x_min,       x_min + (M-1.0)*h, M)
x_B = np.linspace(x_min+0.5*h, x_min + (M-0.5)*h, M)

# %%
plt.plot(x_E, E_arr[700,:])
plt.plot(x_B, B_arr[700,:])
plt.show()

# %%
for n in range(0,N,100):
    plt.figure(figsize=(8,6))
    plt.title(f"n = {n:05d}/N")
    plt.plot(x_E, E_arr[n,:], label="E_y")
    plt.plot(x_B, B_arr[n,:], label="B_z")
    plt.xlim(x_min, x_max)
    plt.ylim(-2,2)
    plt.xlabel("x")
    plt.ylabel("amplitude")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"animation/{n:05d}.png")
    plt.close()

# %%
