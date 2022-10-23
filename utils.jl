function get_length_atomic(l)
    return l / L0
end

function get_length_metric(l)
    return l * L0
end

function get_a(l=L, n=N)
    return l / (2 * n)
end

function find_row_col(index, size)
    row = (index - 1) ÷ size + 1 # integer division
    col = (index - 1) % size + 1

    return row, col
end

function find_index(row, col, size)
    return (row - 1) * size + col
end

function find_row_col_node(index, n=N)
    return find_row_col(index, 2 * n + 1)
end

function find_row_col_elements(index, n=N)
    return find_row_col(index, 2 * n)
end

function find_index_node(row, col, n=N)
    return find_index(row, col, 2 * n + 1)
end

function find_index_elements(row, col, n=N)
    return find_index(row, col, 2 * n)
end

function get_coordinates_nodes(index, l=L, n=N)
    a = get_a(l, n)
    row, col = find_row_col_node(index, n)
    return (col - n - 1) * a, (row - n - 1) * a
end

function get_local_index_node(index, element, n=N)
    row, col = find_row_col_node(index, n)
    row_element, col_element = find_row_col_elements(element, n)
    x = row - row_element
    y = col - col_element
    if x > 1 || y > 1 || x < 0 || y < 0
        error("Element does not contain a node")
    end
    return 2 * x + y + 1
end

function get_global_index_node(local_index, element, n=N)
    row_element, col_element = find_row_col_elements(element, n)
    x = (local_index - 1) ÷ 2
    y = (local_index - 1) % 2
    return find_index_node(row_element + x, col_element + y, n)
end

function print_node(index, element, l=L)
    println("Global index: $index")
    println("Local index: $(get_local_index_node(index, element))")
    println("Coordinates: $(get_coordinates_nodes(index, l))")

end

function f1(ξ)
    return 0.5 * (1 - ξ)
end

function f2(ξ)
    return 0.5 * (1 + ξ)
end

function g1(ξ1, ξ2)
    return f1(ξ1) * f1(ξ2)
end

function g2(ξ1, ξ2)
    return f2(ξ1) * f1(ξ2)
end

function g3(ξ1, ξ2)
    return f1(ξ1) * f2(ξ2)
end

function g4(ξ1, ξ2)
    return f2(ξ1) * f2(ξ2)
end

function g(i, ξ1, ξ2)
    return [g1, g2, g3, g4][i](ξ1, ξ2)
end

function get_x_real(ξ1, k, n=N)
    global_index_1 = get_global_index_node(1, k, n)
    global_index_2 = get_global_index_node(2, k, n)
    x1, _ = get_coordinates_nodes(global_index_1)
    x2, _ = get_coordinates_nodes(global_index_2)
    return x1 * f1(ξ1) + x2 * f2(ξ1)
end

function get_y_real(ξ2, k, l=L, n=N)
    global_index_1 = get_global_index_node(3, k, n)
    global_index_3 = get_global_index_node(4, k, n)
    _, y1 = get_coordinates_nodes(global_index_1, l, n)
    _, y3 = get_coordinates_nodes(global_index_3, l, n)
    return y1 * f1(ξ2) + y3 * f2(ξ2)
end

function get_real_coordinates(ξ1, ξ2, k, l=L, n=N)
    return get_x_real(ξ1, k, n), get_y_real(ξ2, k, l, n)
end

