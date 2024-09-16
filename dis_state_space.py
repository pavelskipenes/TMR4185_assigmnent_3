import numpy as np
from scipy.linalg import expm


def dis_state_space(
        mass: float,
        stiffness: float,
        damping: float,
        forces: list[float],
        time_step: float,
        displacement_initial: float,
        velocity_initial: float):

    A = np.array([[0, 1], [-stiffness/mass, -damping/mass]])
    B = np.array([[0], [1]])

    A_d = expm(A * time_step)
    B_d = np.linalg.inv(A) @ (A_d - np.eye(A.shape[0])) @ B

    states = np.array([[displacement_initial], [velocity_initial]])

    displacements = []

    for force in forces:
        states = A_d @ states + B_d * force
        displacements.append(states[0, 0])

    return displacements
