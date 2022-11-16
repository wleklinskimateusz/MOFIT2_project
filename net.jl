function find_row_col(index::Int64, size::Int64)::Tuple{Int64,Int64}
    row = (index - 1) ÷ size + 1 # integer division
    col = (index - 1) % size + 1

    return row, col
end

function find_index(row::Int64, col::Int64, size::Int64)::Int64
    return (row - 1) * size + col
end

"Return row and column of a node given by an index."
function find_row_col_node(index::Int64, n::Int16)::Tuple{Int64,Int64}
    return find_row_col(index, 2 * n + 1)
end

"Return row and column of an element given by an index."
function find_row_col_elements(index::Int64, n::Int16)::Tuple{Int64,Int64}
    return find_row_col(index, 2 * n)
end

"Return index of a nodebased on it's row and column."
function find_index_node(row::Int64, col::Int64, n::Int16)::Int64
    return find_index(row, col, 2 * n + 1)
end

"Return index of an element based on it's row and column."
function find_index_elements(row::Int64, col::Int64, n::Int16)
    return find_index(row, col, 2 * n)
end

"Return coordinates of a node given it's index."
function get_coordinates_nodes(index::Int64, l::Float64, n::Int16)
    row::Int64, col::Int64 = find_row_col_node(index, n)
    a = get_a(l, n)
    return (col - n - 1) * a, (row - n - 1) * a
end

"Find local node in an element given it's global index."
function get_local_index_node(index::Int64, element::Int64, n::Int16)::Int64
    row, col = find_row_col_node(index, n)
    row_element, col_element = find_row_col_elements(element, n)
    x = row - row_element
    y = col - col_element
    if x > 1 || y > 1 || x < 0 || y < 0
        error("Element does not contain a node")
    end
    return 2 * x + y + 1
end

"Find global node given it's local index and element number."
function get_global_index_node(local_index::Int64, element::Int64, n::Int16)::Int64
    row_element, col_element = find_row_col_elements(element, n)
    x = (local_index - 1) ÷ 2
    y = (local_index - 1) % 2
    return find_index_node(row_element + x, col_element + y, n)
end

"Compute 'x' value in real space based on coordinates in x direction in frame of reference."
function get_x_real(ξ1::Float64, k::Int64, l::Float64, n::Int16)::Float64
    x1, _ = get_coordinates_nodes(get_global_index_node(1, k, n), l, n)
    x2, _ = get_coordinates_nodes(get_global_index_node(2, k, n), l, n)
    return x1 * f1(ξ1) + x2 * f2(ξ1)
end

"Compute 'y' value in real space based on coordinates in y direction in frame of reference."
function get_y_real(ξ2::Float64, k::Int64, l::Float64, n::Int16)::Float64
    _, y1 = get_coordinates_nodes(get_global_index_node(1, k, n), l, n)
    _, y3 = get_coordinates_nodes(get_global_index_node(3, k, n), l, n)
    return y1 * f1(ξ2) + y3 * f2(ξ2)
end

"Compute 'x, y' value in real space based on coordinates in x,y directions in frame of reference."
function get_real_coordinates(ξ1::Float64, ξ2::Float64, k::Int64, l::Float64, n::Int16)::Tuple{Float64,Float64}
    return get_x_real(ξ1, k, l, n), get_y_real(ξ2, k, l, n)
end


"Print node's global index, local nidex and coordinates given it's global index."
function print_node(index::Int64, element::Int64, n::Int16, l::Float64)::Nothing
    println("Global index: $index")
    println("Local index: $(get_local_index_node(index, element, n))")
    println("Coordinates: $(get_coordinates_nodes(index, l, n))")
end

"Print information about all nodes in an element."
function print_element(k::Int64, n::Int16, l::Float64)::Nothing
    println("element $k")
    for local_node in 1:4
        global_node = get_global_index_node(local_node, k, n)
        println("local node: $local_node, global node: $global_node")
        print_node(global_node, k, n, l)
        println()
    end
end

"Print information about all elements."
function print_elements(n::Int16, l::Float64)::Nothing
    for element_idx in 1:4*n^2
        print_element(element_idx, n, l)
        println()
    end
end
