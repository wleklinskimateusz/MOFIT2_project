from const import L, N, M, OMEGA
from utils import find_row_col_node, get_coordinates_nodes, global_index_node, print_node, g
import numpy as np
import matplotlib.pyplot as plt


def print_elements():
    for element_idx in range(4*N**2):
        print(f"[ELEMENT]: {element_idx}")
        for local_node in range(4):
            global_node = global_index_node(local_node, element_idx)
            print_node(global_node, element_idx)
            print()


def get_psi_nodes(x, y, m=M, omega=OMEGA):
    return np.exp(-m*omega/2*(x**2 + y**2))


def get_psi_teo():
    psi = np.zeros((20 * 2*N + 1, 20 * 2*N + 1))
    for i in range(psi.shape[0]):
        for j in range(psi.shape[1]):
            psi[i, j] = get_psi_nodes(
                i * 0.1 * L / (4 * N) - L/2, j * 0.1 * L / (4 * N) - L/2)
    return psi


def get_psi(k: int, psi: np.ndarray, ksi_1: float, ksi_2: float) -> np.ndarray:
    output = 0
    for i in range(4):
        idx = global_index_node(i, k)
        row, col = find_row_col_node(idx)
        local_psi = psi[row*20, col*20]
        output += local_psi * g(i, ksi_1, ksi_2)
    return output


def populate_nodes(psi):
    for i in range((2*N+1)**2):
        x, y = get_coordinates_nodes(i)
        row, col = find_row_col_node(i)
        psi[row*20, col*20] = get_psi_nodes(x, y)


def calculate_psi_element(element: int, psi: np.ndarray):
    ksi_1 = np.arange(-1, 1, 0.1)
    ksi_2 = np.arange(-1, 1, 0.1)
    idx = global_index_node(0, element)
    row, col = find_row_col_node(idx)
    for i, k1 in enumerate(ksi_1):
        for j, k2 in enumerate(ksi_2):
            if ((i, j) in [(0, 0), (0, 20), (20, 0), (20, 20)]):
                continue
            psi[20*row+j, 20*col+i] = get_psi(element, psi, k1, k2)


def plot_psi(psi: np.ndarray):
    plt.imshow(psi)
    plt.colorbar()
    plt.show()


def main():
    print_elements()  # ex1

    psi = np.zeros((20 * 2*N + 1, 20 * 2*N + 1))
    populate_nodes(psi)

    for element in range(4*N**2):
        calculate_psi_element(element, psi)

    plot_psi(psi)

    plot_psi(get_psi_teo())


if __name__ == "__main__":
    main()
