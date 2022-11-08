function get_length_atomic(l)
    return l / L0
end

function get_length_metric(l)
    return l * L0
end

function get_a(l=L)
    return l / (2 * N)
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
    row, col = find_row_col_node(index, n)
    return (col - n - 1) * get_a(l), (row - n - 1) * get_a(l)
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

function get_x_real(ξ1, k, l=L, n=N)
    global_index_1 = get_global_index_node(1, k, n)
    global_index_2 = get_global_index_node(2, k, n)
    x1, _ = get_coordinates_nodes(global_index_1, l, n)
    x2, _ = get_coordinates_nodes(global_index_2, l, n)
    return x1 * f1(ξ1) + x2 * f2(ξ1)
end

function get_y_real(ξ2, k, l=L, n=N)
    global_index_1 = get_global_index_node(1, k, n)
    global_index_3 = get_global_index_node(3, k, n)
    _, y1 = get_coordinates_nodes(global_index_1, l, n)
    _, y3 = get_coordinates_nodes(global_index_3, l, n)
    return y1 * f1(ξ2) + y3 * f2(ξ2)
end

function get_real_coordinates(ξ1, ξ2, k, l=L, n=N)
    return get_x_real(ξ1, k, l, n), get_y_real(ξ2, k, l, n)
end


function get_s_element(i, j, len=L)
    output = 0
    for l in 1:3
        for n in 1:3
            output += w[l] * w[n] * g(j, p[l], p[n]) * g(i, p[l], p[n])
        end
    end
    return output * get_a(len)^2 / 4
end


function get_matrix(get_matrix_element)
    output = zeros(4, 4)
    for i in 1:4
        for j in 1:4
            output[i, j] = get_matrix_element(i, j)
        end
    end
    return output
end

function dev_g_ξ_1(j, x1, x2)
    return (g(j, x1, x2 + ε) - g(j, x1, x2 - ε)) / (2 * ε)
end

function dev_g_ξ_2(j, x1, x2)
    return (g(j, x1 + ε, x2) - g(j, x1 - ε, x2)) / (2 * ε)
end

function get_t_element(i, j)
    output = 0
    for l in 1:3
        for n in 1:3
            output += w[l] * w[n] * (dev_g_ξ_1(j, p[l], p[n]) * dev_g_ξ_1(i, p[l], p[n]) + dev_g_ξ_2(j, p[l], p[n]) * dev_g_ξ_2(i, p[l], p[n]))
        end
    end
    return output / (2 * M)
end

function get_v_element(i, j, k, len=L)
    output = 0
    for l in 1:3
        for n in 1:3
            output += w[l] * w[n] * g(j, p[l], p[n]) * g(i, p[l], p[n]) * (get_x_real(p[l], k, len)^2 + get_y_real(p[n], k, len)^2)
        end
    end
    return get_a(len)^2 / 4 * M * OMEGA^2 / 2 * output

end

function get_global_matrix(l=L)
    H = zeros((2 * N + 1)^2, (2 * N + 1)^2)
    S = zeros((2 * N + 1)^2, (2 * N + 1)^2)
    for k in 1:4*N^2
        for i in 1:4
            for j in 1:4
                global_index_1 = get_global_index_node(i, k)
                global_index_2 = get_global_index_node(j, k)
                H[global_index_1, global_index_2] += get_v_element(i, j, k, l) + get_t_element(i, j)
                S[global_index_1, global_index_2] += get_s_element(i, j, l)
            end
        end
    end
    return H, S
end

function handle_removing_things(H, S, i)
    for j in 1:(2*N+1)^2
        H[i, j] = 0
        H[j, i] = 0
        S[i, j] = 0
        S[j, i] = 0
    end
    S[i, i] = 1
    H[i, i] = -bitwaPodGrunwaldem
end

function get_edge_indexes()
    output = []
    i = 1
    while i <= 2 * N + 1
        push!(output, i)
        i += 1
    end
    while i <= (2 * N + 1)^2 - (2 * N + 1)
        push!(output, i)
        i += 2 * N
        push!(output, i)
        i += 1
    end
    while i <= (2 * N + 1)^2
        push!(output, i)
        i += 1
    end
    return output
end

function print_table(tab, n)
    println("New TABLE:")
    for i in 1:n
        for j in 1:n
            print(tab[i, j], "\t")
        end
        println()
    end
end

function save_matrix_to_file(A, size, filename)
    open("output/$filename", "w") do io
        for i in 1:size
            for j in 1:size
                write(io, "$(A[i, j])\t")
            end
            write(io, "\n")
        end
    end
end

