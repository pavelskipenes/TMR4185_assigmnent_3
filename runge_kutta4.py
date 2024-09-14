import numpy as np

def runge_kutta4(
    mass: float,
    spring: float,
    damping: float,
    forces: list[float],
    time_step: float,
    displacement_initial: float,
    velocity_initial: float
) -> list[float]:

    def get_states_dot(state: np.ndarray, input: float, damper: float, spring: float, mass: float) -> np.ndarray:
        return np.array([state[1], (input - spring * state[0] - damper * state[1]) / mass])

    def get_next_state(time_step: float, state: np.ndarray, input: float, damper: float, spring: float, mass: float) -> np.ndarray:
        K0 = get_states_dot(state, input, damper, spring, mass)
        K1 = get_states_dot(state + (time_step / 2) * K0, input, damper, spring, mass)
        K2 = get_states_dot(state + (time_step / 2) * K1, input, damper, spring, mass)
        K3 = get_states_dot(state + time_step * K2, input, damper, spring, mass)
        return state + (time_step / 6) * (K0 + 2 * K1 + 2 * K2 + K3)

    state = np.array([displacement_initial, velocity_initial])
    
    states = np.array([get_next_state(time_step, state, force, damping, spring, mass) for force in forces])
    
    return states[:,0]