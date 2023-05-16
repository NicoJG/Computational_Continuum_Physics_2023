import numpy as np
import matplotlib.pyplot as plt

# For latex interpretation of the figures
#plt.rcParams.update({
#    "text.usetex": True,
#    "font.family": "Computer Modern",
#    "font.size": 11.0, # 11pt is fontsize of captions in the report
#})

# Figure setup
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(1,1,1)

# Iterate over data files
files = ["saved_data/psi_time_and_frequency_harmonic.csv"]
for i, filename in enumerate(files):

    # Extract data
    data = np.genfromtxt(filename, unpack=True, delimiter=",", skip_header=1)
    time, omega, abs_psi_t, arg_psi_t, abs_psi_omega, arg_psi_omega = data
    i_sorted = np.argsort(omega)

    # find the peak:
    idx_peak = np.argmax(abs_psi_omega)
    omega_peak = omega[idx_peak]
    print(f"Omega of the highest peak: {omega_peak:.5f}")

    # Add spectrum to plot
    ax.plot(omega[i_sorted], abs_psi_omega[i_sorted]/len(omega),
            ["b", "r"][i], lw=2,
            label="$t_\mathrm{max}=" + f"{round(time[-1])}$")

# Final figure formatting
ax.set_xlabel("$\omega$")
ax.set_ylabel("$\psi(0.1, 0, \omega)$")
ax.set_xlim(-10, 2)
#ax.legend()
ax.grid(which="both", axis="x")
ax.grid(axis="y")
ax.minorticks_on()
fig.tight_layout()
fig.savefig("figures/spectrum.png")
fig.savefig("figures/spectrum.pdf")

