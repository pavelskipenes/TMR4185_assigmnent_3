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

