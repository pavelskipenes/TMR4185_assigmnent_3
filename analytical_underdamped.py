import numpy as np

def voldsfunksjon(u0, udot0, k, P0, omega, omega_d, omega_0, t, epsilon, beta, zeta):
    # Løser lign i rekkefølge slik at verdier finnes
    stor_u = (P0/k)/np.sqrt((1-beta**2)**2 + (2*zeta*beta)**2)
    phi = np.arctan((-2*zeta*beta)/(1-beta**2))
    theta = np.arctan((udot0 + zeta*omega_0*(u0 - stor_u*np.cos(epsilon + phi)) + omega*stor_u*np.sin(epsilon + phi))
                      /(omega_d*(u0 - stor_u*np.cos(epsilon + phi))))
    R = ((u0 - stor_u*np.cos(epsilon + phi))
         /np.cos(theta))
    #løse oppgave
    u = ((np.exp(-zeta*omega_0*t)*R*np.cos(omega_d*t - theta))
         +(stor_u*np.cos(omega*t + epsilon + phi)))
    return u