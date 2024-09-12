import numpy as np
import matplotlib.pyplot as plt

from const_avg_acc import const_avg_acc
from runge_kutta4 import runge_kutta4
from dis_state_space import dis_state_space
from analytical_underdamped import voldsfunksjon

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
uas = voldsfunksjon(c, m, u0, udot0, k, P[0], omega, t, eps1)

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
plt.plot(t, uas, 'y--', label='VOLD')
plt.legend()
plt.grid(True)
plt.xlabel('time, s')
plt.ylabel('u, m')
# Dette burde ha dekt deloppgave 1 innen kode.

#går så over til del 2