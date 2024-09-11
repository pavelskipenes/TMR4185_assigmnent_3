import numpy as np
import matplotlib.pyplot as plt

from const_avg_acc import const_avg_acc
from runge_kutta4 import runge_kutta4
from dis_state_space import dis_state_space

# Input

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

# Student code called here!
u = const_avg_acc(m, k, c, P, h, u0, udot0)
urk = runge_kutta4(m, k, c, P, h, u0, udot0)
uss = dis_state_space(m, k, c, P, h, u0, udot0)

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
plt.legend()
plt.grid(True)
plt.xlabel('time, s')
plt.ylabel('u, m')
