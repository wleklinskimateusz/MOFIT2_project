from const import L, L0, N


def get_length_atomic(l: float) -> float:
    return l / L0


def get_length_metric(l: float) -> float:
    return l * L0


def get_a(l: float = L, n: int = N):
    return l / (2 * n)


def find_row_col(index: int, size: int) -> tuple[int]:
    row = index // size
    col = index % size
    return row, col


def find_index(row: int, col: int, size: int) -> int:
    return row * size + col


def find_row_col_node(index: int, n: int = N) -> tuple[int]:
    return find_row_col(index, 2 * n + 1)


def find_row_col_elements(index: int, n: int = N) -> tuple[int]:
    return find_row_col(index, 2 * n)


def find_index_node(row: int, col: int, n: int = N) -> int:
    return find_index(row, col, 2 * n + 1)


def find_index_elements(row: int, col: int, n: int = N) -> int:
    return find_index(row, col, 2 * n)


def get_coordinates_nodes(index: int, l: float = L, n: int = N) -> tuple[float]:
    row, col = find_row_col_node(index, n)
    a = get_a(l, n)
    x = (col - n) * a
    y = (row - n) * a
    return x, y


def local_index_node(index: int,  element: int, n: int = N) -> int:
    row, col = find_row_col_node(index, n)
    row_element, col_element = find_row_col_elements(element, n)
    x = row - row_element
    y = col - col_element
    if x > 1 or y > 1:
        raise Exception("element does not contain node")
    return 2 * x + y  # 1, 2, 3, 4


def global_index_node(local_index: int, element: int, n: int = N) -> int:
    row_element, col_element = find_row_col_elements(element, n)
    row = row_element + local_index // 2
    col = col_element + local_index % 2
    return find_index_node(row, col, n)


def print_node(index: int,  element: int, l: float = L):
    print(f"global index: {index}")
    print(f"local index: {local_index_node(index, element)}")
    print(f"coordinates: {get_coordinates_nodes(index, l)} [nm]")


def f1(ksi: float) -> float:
    return 0.5 * (1 - ksi)


def f2(ksi: float) -> float:
    return 0.5 * (1 + ksi)


def g1(ksi_1: float, ksi_2: float) -> float:
    return f1(ksi_1) * f1(ksi_2)


def g2(ksi_1: float, ksi_2: float) -> float:
    return f2(ksi_1) * f1(ksi_2)


def g3(ksi_1: float, ksi_2: float) -> float:
    return f1(ksi_1) * f2(ksi_2)


def g4(ksi_1: float, ksi_2: float) -> float:
    return f2(ksi_1) * f2(ksi_2)


def g(i: int, ksi_1: float, ksi_2: float) -> float:
    return [g1, g2, g3, g4][i](ksi_1, ksi_2)


def get_x_real(ksi_1: float, k: int, n: int = N):
    global_index_1 = global_index_node(1, k, n)
    global_index_2 = global_index_node(2, k, n)
    x1, _ = get_coordinates_nodes(global_index_1, L, n)
    x2, _ = get_coordinates_nodes(global_index_2, L, n)
    return x1*f1(ksi_1) + x2*f2(ksi_1)


def get_y_real(ksi_2: float, k: int, n: int = N):
    global_index_1 = global_index_node(1, k, n)
    global_index_3 = global_index_node(3, k, n)
    _, y1 = get_coordinates_nodes(global_index_1, L, n)
    _, y3 = get_coordinates_nodes(global_index_3, L, n)
    return y1*f1(ksi_2) + y3*f2(ksi_2)


def get_real_coordinates(ksi_1: float, ksi_2: float, k: int, n: int = N):
    return get_x_real(ksi_1, k, n), get_y_real(ksi_2, k, n)
