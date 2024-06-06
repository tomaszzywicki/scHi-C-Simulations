import numpy as np

def lj_potential(distance, epsilon=1, sigma=1):
    """
    Oblicza potencjał Lennarda-Jonesa dla danej odległości między dwoma cząsteczkami.
    
    :param distance: Odległość między cząsteczkami.
    :param epsilon: Parametr epsilon potencjału (domyślnie 1).
    :param sigma: Parametr sigma potencjału (domyślnie 1).
    :return: Wartość potencjału Lennarda-Jonesa.
    """
    if distance == 0:
        return 0  # Jeśli odległość wynosi zero, potencjał również wynosi zero
    else:
        r6 = (sigma / distance) ** 6
        return 4 * epsilon * (r6 ** 2 - r6) if r6 != 0 else 0

def lj_energy(walk, epsilon=1, sigma=1):
    """
    Oblicza całkowitą energię potencjalną spaceru na podstawie potencjału Lennarda-Jonesa.
    
    :param walk: Lista punktów spaceru.
    :param epsilon: Parametr epsilon potencjału (domyślnie 1).
    :param sigma: Parametr sigma potencjału (domyślnie 1).
    :return: Całkowita energia potencjalna spaceru.
    """
    total_energy = 0
    n = len(walk)
    for i in range(n):
        for j in range(i+1, n):
            distance = np.linalg.norm(np.array(walk[i]) - np.array(walk[j]))
            total_energy += lj_potential(distance, epsilon, sigma)
    return total_energy