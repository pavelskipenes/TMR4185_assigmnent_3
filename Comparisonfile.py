import numpy as np
import matplotlib.pyplot as plt

from const_avg_acc import const_avg_acc
from runge_kutta4 import runge_kutta4
from dis_state_space import dis_state_space
from analytical_underdamped import analytical_underdamped

# Input dynamic system characteristics
m = 1  # mass, kg
k = 1  # stiffness, N/m
c = 0.5  # damping, N/s
u0 = 0.2  # initial position, m
udot0 = 0.1  # initial velocity, m/s
loadtype = 2  # see if statements below for definition

h = 0.1  # time increment, s
tmax = 400  # time duration, s
t = (np.arange(0, tmax, h))
P = np.zeros(len(t))  # initialize load vector to zero

#Verdier spess for harmonic load
omega = 0.4  # frequency
amp = 1  # amplitude
eps1 = np.pi  # phase angle
P = amp * np.cos(omega * t + eps1)

# Student code called here!
u = const_avg_acc(m, k, c, P, h, u0, udot0)
urk = runge_kutta4(m, k, c, P, h, u0, udot0)
uss = dis_state_space(m, k, c, P, h, u0, udot0)
uas = analytical_underdamped(c, m, u0, udot0, k, P[0], omega, t, eps1)

# Plot
plt.figure(1)
plt.subplot(211)
plt.plot(t, P)
plt.grid(True)
plt.xlabel('time, s')
plt.ylabel('P, N')
plt.xlim((0, tmax))

plt.subplot(212)
plt.plot(t, u, 'k', label='CAA')
plt.plot(t, urk, 'b--', label='RK')
plt.plot(t, uss, 'r--', label='DSS')
plt.plot(t, uas, 'y--', label='AU')
plt.legend()
plt.grid(True)
plt.xlabel('time, s')
plt.ylabel('u, m')
# Dette burde ha dekt deloppgave 1 innen kode.

#går så over til del 2
# Parameters for different cases
cases = {
    "Resonance": {"c": 0.5, "m": 1.0, "u0": 0.0, "udot0": 0.0, "k": 1.0, "P0": 1.0, "omega": 1.0, "epsilon": 0.0},
    "Stiffness-Dominated": {"c": 0.5, "m": 1.0, "u0": 0.0, "udot0": 0.0, "k": 10.0, "P0": 1.0, "omega": 0.5, "epsilon": 0.0},
    "Inertia-Dominated": {"c": 0.5, "m": 1.0, "u0": 0.0, "udot0": 0.0, "k": 1.0, "P0": 1.0, "omega": 2.0, "epsilon": 0.0}
}

# Time vector
t = np.linspace(0, 20, 1000)

# Storage for analysis
steady_state_amplitudes = {}
transient_durations = {}

plt.figure(figsize=(14, 10))

for case, params in cases.items():
    u = analytical_underdamped(**params, t=t)
    plt.plot(t, u, label=f"{case}")

    # Calculate steady-state amplitude
    steady_state_amplitude = np.max(u[-100:])  # Last 100 points as steady-state
    steady_state_amplitudes[case] = steady_state_amplitude

    # Find transient duration
    threshold = 0.05 * steady_state_amplitude
    transient_duration = np.max(t[u > threshold]) if np.any(u > threshold) else np.nan
    transient_durations[case] = transient_duration

plt.title("System Response for Different Types of Systems")
plt.xlabel("Time")
plt.ylabel("Displacement")
plt.legend()
plt.grid(True)
plt.show()

# Print the results
print("Steady-State Amplitudes:")
for case, amplitude in steady_state_amplitudes.items():
    print(f"{case}: {amplitude:.2f}")

print("\nTransient Durations (time when transient response < 5% of steady-state amplitude):")
for case, duration in transient_durations.items():
    print(f"{case}: {duration:.2f} seconds")
