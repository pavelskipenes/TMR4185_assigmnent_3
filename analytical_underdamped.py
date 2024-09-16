import numpy as np


def analytical_underdamped(c, m, u0, udot0, k, P0, omega, time, epsilon=0):
    omega_0 = np.sqrt(k/m)
    zeta = c/(2*np.sqrt(m*k))
    beta = omega/omega_0
    omega_d = omega_0*np.sqrt(1-zeta**2)

    U = (P0/k)/np.sqrt((1-beta**2)**2 + (2*zeta*beta)**2)
    phi = np.arctan((-2*zeta*beta)/(1-beta**2))

    theta = np.arctan((udot0 + zeta*omega_0*(u0 - U*np.cos(epsilon + phi)) + omega*U*np.sin(epsilon + phi))
                      / (omega_d*(u0 - U*np.cos(epsilon + phi))))

    R = (u0 - U*np.cos(epsilon + phi)) / np.cos(theta)

    displacements = []
    for t in time:
        displacement_transient = R * np.exp(-zeta*omega_0*t)*np.cos(omega_d*t - theta)
        displacement_ss = U*np.cos(omega*t + epsilon + phi)
        displacements.append(displacement_ss + displacement_transient)
    return displacements
