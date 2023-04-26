import numpy as np
import matplotlib.pyplot as plt

# For latex interpretation of the figures
#plt.rcParams.update({
#    "text.usetex": True,
#    "font.family": "Computer Modern",
#    "font.size": 11.0, # 11pt is fontsize of captions in the report
#})

# time, omega, Re(psi(t)), Im(psi(t)), Re(psi(omega)), Im(psi(omega))
data = np.genfromtxt("data/psi_time_and_frequency.csv", unpack=True, delimiter=",", skip_header=1)
time, omega, abs_psi_t, arg_psi_t, abs_psi_omega, arg_psi_omega = data

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(omega, abs_psi_omega)
ax.set_xlabel("$\omega$")
ax.set_ylabel("$\psi(0, 0.1, \omega)$ (a.u.)")
ax.set_xlim(-1, 5)
ax.grid()
fig.savefig("figures/spectrum.png")

# find the peak:
idx_peak = np.argmax(abs_psi_omega)
omega_peak = omega[idx_peak]
print(f"Omega of the highest peak: {omega_peak:.5f}")