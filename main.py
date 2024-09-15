import numpy as np
import matplotlib.pyplot as plt

from const_avg_acc import const_avg_acc
from runge_kutta4 import runge_kutta4
from dis_state_space import dis_state_space
from analytical_underdamped import analytical_underdamped_split
from analytical_underdamped import analytical_underdamped

# Input

# Input dynamic system characteristics
m = 1  # mass, kg
k = 1  # stiffness, N/m
c = 0.5  # damping, N/s
u0 = 0.2  # initial position, m
udot0 = 0.1  # initial velocity, m/s
loadtype = 2  # see if statements below for definition

h = 0.1  # time increment, s
tmax = 25  # time duration, s
t = (np.arange(0, tmax, h))
P = np.zeros(len(t))  # initialize load vector to zero

# Formulate load vector
if (loadtype == 0):  # harmonic load
    omega = 0.4  # frequency
    amp = 1  # amplitude
    eps1 = np.pi  # phase angle
    P = amp * np.cos(omega * t + eps1)  # load vector

elif (loadtype == 1):  # general periodic load
    ncomponents = 3  # number of load components
    omegas = np.array([0.9, 1, 1.1, 0.25, 1/6])  # load frequencies
    amps = [1, 1, 1, 1, 4]  # load amplitudes
    epss = np.array([0, np.pi/4, 0, np.pi, 3*np.pi/4])

    for ii in range(0, ncomponents):
        P = P + amps[ii] * np.cos(omegas[ii] * t + epss[ii])

else:  # (loadtype == 3) short duration
    Pmax = 10  # maximum load, N
    t1 = 1  # time duration of load, s
    omega = 2 * np.pi / (2 * t1)  # frequency for sine shape
    indt1 = sum(t <= t1)  # find end of the impulse
    P[0:indt1] = Pmax * np.sin(omega * t[0:indt1])

u = const_avg_acc(m, k, c, P, h, u0, udot0)
urk = runge_kutta4(m, k, c, P, h, u0, udot0)
uss = dis_state_space(m, k, c, P, h, u0, udot0)
analytical_stationary, analytical_transient = analytical_underdamped_split(c, m, u0,udot0, k,10,omega, t, 0)
analytical = analytical_underdamped(c, m, u0,udot0, k,10,omega, t, 0)

# ignore plotting the impulse load as it's not that interesting
#plt.figure(1)
#plt.subplot(211)
#plt.plot(t, P)
#plt.grid(True)
#plt.xlabel('time, s')
#plt.ylabel('P, N')
#plt.xlim((0, tmax))

#plt.subplot(212)
plt.plot(t, u, 'k', label='CAA')
plt.plot(t, urk, 'b--', label='RK')
plt.plot(t, uss, 'r--', label='DSS')
plt.plot(t, analytical_transient, 'g-', label='analytical transient')
plt.plot(t, analytical_stationary, 'y-', label='analytical stationary')
plt.plot(t, analytical, 'y-', label='analytical')
plt.legend()
plt.grid(True)
plt.xlabel('time, s')
plt.ylabel('u, m')
plt.savefig("plots/analytical_and_numerical_aproximations.png")
plt.clf()

plt.plot(t, analytical_transient, 'g-', label='analytical transient')
plt.plot(t, analytical_stationary, 'y-', label='analytical stationary')
plt.legend()
plt.grid(True)
plt.xlabel('time, s')
plt.ylabel('u, m')
plt.savefig("plots/analytical_analysis.png")
plt.clf()

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
    u = analytical_underdamped(**params, time=t)
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
plt.savefig("plots/system_response_for_different_types_of_systems.png")
plt.clf()

# Print the results
print("Steady-State Amplitudes:")
for case, amplitude in steady_state_amplitudes.items():
    print(f"{case}: {amplitude:.2f}")

print("\nTransient Durations (time when transient response < 5% of steady-state amplitude):")
for case, duration in transient_durations.items():
    print(f"{case}: {duration:.2f} seconds")