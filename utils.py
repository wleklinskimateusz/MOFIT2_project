def get_a(l, n):
    return l / (2 * n)


def find_row_col(index: int, size: int) -> tuple[int]:
    row = index // size
    col = index % size
    return row, col


def find_index(row: int, col: int, size: int) -> int:
    return row * size + col


def find_row_col_node(index: int, n: int) -> tuple[int]:
    return find_row_col(index, 2 * n + 1)


def find_row_col_elements(index: int, n: int) -> tuple[int]:
    return find_row_col(index, 2 * n)


def find_index_node(row: int, col: int, n: int) -> int:
    return find_index(row, col, 2 * n + 1)


def find_index_elements(row: int, col: int, n: int) -> int:
    return find_index(row, col, 2 * n)


def get_coordinates_nodes(index: int, n: int, l: float) -> tuple[float]:
    row, col = find_row_col_node(index, n)
    a = get_a(l, n)
    x = (col - n) * a
    y = (row - n) * a
    return x, y


def local_index_node(index: int, n: int, element: int) -> int:
    row, col = find_row_col_node(index, n)
    row_element, col_element = find_row_col_elements(element, n)
    x = row - row_element
    y = col - col_element
    if x > 1 or y > 1:
        raise Exception("element does not contain node")
    return 2 * x + y


def global_index_node(local_index: int, n: int, element: int) -> int:
    row_element, col_element = find_row_col_elements(element, n)
    row = row_element + local_index // 2
    col = col_element + local_index % 2
    return find_index_node(row, col, n)


def print_node(index: int, n: int, l: float, element: int):
    print(f"global index: {index}")
    print(f"local index: {local_index_node(index, n, element)}")
    print(f"coordinates: {get_coordinates_nodes(index, n, l)} [nm]")
