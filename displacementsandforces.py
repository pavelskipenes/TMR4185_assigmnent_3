import numpy as np
import matplotlib.pyplot as plt

# Time settings
t = np.linspace(0, 30, 1000)
omega_0 = 1  # Natural frequency (rad/s)

# External force (resonance condition: omega ≈ omega_0)
P = np.cos(omega_0 * t)

# System parameters
zeta_values = [0.01, 0.1, 1]  # Different damping ratios
methods = ['CAA', 'RK', 'DSS']  # Numerical methods

# Function for displacement response (mock-up)
def displacement_response(t, omega_0, zeta):
    return np.exp(-zeta * t) * np.sin(omega_0 * t)

# Prepare figure and axes
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Plot the external force
ax1.plot(t, P, color='blue', label='t1 = 1')
ax1.set_title(r'System at resonance, $\omega \approx \omega_0$')
ax1.set_xlabel('time, s')
ax1.set_ylabel(r'$P$, N')
ax1.legend()
ax1.grid(True)

# Plot the displacement responses for different methods and zeta values
colors = ['cyan', 'magenta', 'yellow']
linestyles = ['--', '-.', ':']
for i, zeta in enumerate(zeta_values):
    for j, method in enumerate(methods):
        u = displacement_response(t, omega_0, zeta)
        ax2.plot(t, u, linestyle=linestyles[j], color=colors[i],
                 label=f'{method}, ζ = {zeta}')

ax2.set_xlabel('time, s')
ax2.set_ylabel(r'$u$, m')
ax2.legend()
ax2.grid(True)

# Show plot
plt.tight_layout()
plt.show()
