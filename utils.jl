using Printf

function get_length_atomic(l::Float64)::Float64
    return l / L0
end

function get_length_metric(l::Float64)::Float64
    return l * L0
end

function get_a(l::Float64=L)::Float64
    return l / (2 * N)
end


function find_row_col(index::Int64, size::Int64)::Tuple{Int64,Int64}
    row = (index - 1) ÷ size + 1 # integer division
    col = (index - 1) % size + 1

    return row, col
end

function find_index(row::Int64, col::Int64, size::Int64)::Int64
    return (row - 1) * size + col
end

"Return row and column of a node given by an index."
function find_row_col_node(index::Int64, n::Int16=N)::Tuple{Int64,Int64}
    return find_row_col(index, 2 * n + 1)
end

"Return row and column of an element given by an index."
function find_row_col_elements(index::Int64, n::Int16=N)::Tuple{Int64,Int64}
    return find_row_col(index, 2 * n)
end

"Return index of a nodebased on it's row and column."
function find_index_node(row::Int64, col::Int64, n::Int16=N)::Int64
    return find_index(row, col, 2 * n + 1)
end

"Return index of an element based on it's row and column."
function find_index_elements(row::Int64, col::Int64, n::Int16=N)
    return find_index(row, col, 2 * n)
end

"Return coordinates of a node given it's index."
function get_coordinates_nodes(index, l=L, n=N)
    row::Int64, col::Int64 = find_row_col_node(index, n)
    return (col - n - 1) * get_a(l), (row - n - 1) * get_a(l)
end

"Find local node in an element given it's global index."
function get_local_index_node(index::Int64, element::Int64, n::Int16=N)::Int64
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
function get_global_index_node(local_index::Int64, element::Int64, n::Int16=N)::Int64
    row_element, col_element = find_row_col_elements(element, n)
    x = (local_index - 1) ÷ 2
    y = (local_index - 1) % 2
    return find_index_node(row_element + x, col_element + y, n)
end

"Print node's global index, local nidex and coordinates given it's global index."
function print_node(index::Int64, element::Int64, l::Float64=L)::Nothing
    println("Global index: $index")
    println("Local index: $(get_local_index_node(index, element))")
    println("Coordinates: $(get_coordinates_nodes(index, l))")
end

"Print information about all nodes in an element."
function print_element(k::Int64)::Nothing
    println("element $k")
    for local_node in 1:4
        global_node = get_global_index_node(local_node, k)
        println("local node: $local_node, global node: $global_node")
        print_node(global_node, k)
        println()
    end
end

"Print information about all elements."
function print_elements(n::Int16=N)::Nothing
    for element_idx in 1:4*n^2
        print_element(element_idx)
        println()
    end
end

# --------------------- BASE FUNCTIONS ---------------------
function f1(ξ::Float64)::Float64
    return 0.5 * (1 - ξ)
end

function f2(ξ::Float64)::Float64
    return 0.5 * (1 + ξ)
end

function g1(ξ1::Float64, ξ2::Float64)::Float64
    return f1(ξ1) * f1(ξ2)
end

function g2(ξ1::Float64, ξ2::Float64)::Float64
    return f2(ξ1) * f1(ξ2)
end

function g3(ξ1::Float64, ξ2::Float64)::Float64
    return f1(ξ1) * f2(ξ2)
end

function g4(ξ1::Float64, ξ2::Float64)::Float64
    return f2(ξ1) * f2(ξ2)
end

function g(i::Int64, ξ1::Float64, ξ2::Float64)::Float64
    return [g1, g2, g3, g4][i](ξ1, ξ2)
end

"Compute 'x' value in real space based on coordinates in x direction in frame of reference."
function get_x_real(ξ1::Float64, k::Int64, l::Float64=L, n::Int16=N)::Float64
    x1, _ = get_coordinates_nodes(get_global_index_node(1, k, n), l, n)
    x2, _ = get_coordinates_nodes(get_global_index_node(2, k, n), l, n)
    return x1 * f1(ξ1) + x2 * f2(ξ1)
end

"Compute 'y' value in real space based on coordinates in y direction in frame of reference."
function get_y_real(ξ2::Float64, k::Int64, l::Float64=L, n::Int16=N)::Float64
    _, y1 = get_coordinates_nodes(get_global_index_node(1, k, n), l, n)
    _, y3 = get_coordinates_nodes(get_global_index_node(3, k, n), l, n)
    return y1 * f1(ξ2) + y3 * f2(ξ2)
end

"Compute 'x, y' value in real space based on coordinates in x,y directions in frame of reference."
function get_real_coordinates(ξ1::Float64, ξ2::Float64, k::Int64, l::Float64=L, n::Int16=N)::Tuple{Float64,Float64}
    return get_x_real(ξ1, k, l, n), get_y_real(ξ2, k, l, n)
end

# --------------------- MATRIXES ---------------------
function get_s_element(i::Int64, j::Int64, len::Float64=L)::Float64
    output::Float64 = 0
    for l in 1:3
        for n in 1:3
            output += w[l] * w[n] * g(j, p[l], p[n]) * g(i, p[l], p[n])
        end
    end
    return output * get_a(len)^2 / 4
end

function dev_g_ξ_1(j::Int64, x1::Float64, x2::Float64)::Float64
    return (g(j, x1, x2 + ε) - g(j, x1, x2 - ε)) / (2 * ε)
end

function dev_g_ξ_2(j::Int64, x1::Float64, x2::Float64)::Float64
    return (g(j, x1 + ε, x2) - g(j, x1 - ε, x2)) / (2 * ε)
end

function get_t_element(i::Int64, j::Int64)::Float64
    output = 0
    for l in 1:3
        for n in 1:3
            output += w[l] * w[n] * (dev_g_ξ_1(j, p[l], p[n]) * dev_g_ξ_1(i, p[l], p[n]) + dev_g_ξ_2(j, p[l], p[n]) * dev_g_ξ_2(i, p[l], p[n]))
        end
    end
    return output / (2 * M)
end

"Create local matrix with values given by function 'get_matrix_element'."
function get_matrix(get_matrix_element::Function)::Matrix{Float64}
    output::Matrix{Float64} = zeros(4, 4)
    for i in 1:4
        for j in 1:4
            output[i, j] = get_matrix_element(i, j)
        end
    end
    return output
end

function get_v_element(i::Int64, j::Int64, k::Int64, len::Float64=L)::Float64
    output = 0
    for l in 1:3
        for n in 1:3
            output += w[l] * w[n] * g(j, p[l], p[n]) * g(i, p[l], p[n]) * (get_x_real(p[l], k, len)^2 + get_y_real(p[n], k, len)^2)
        end
    end
    return get_a(len)^2 / 4 * M * OMEGA^2 / 2 * output
end

"Create global S and H matrixes."
function get_global_matrix(len::Float64=L, n::Int16=N)::Tuple{Matrix{Float64},Matrix{Float64}}
    H::Matrix{Float64} = zeros((2 * n + 1)^2, (2 * n + 1)^2)
    S::Matrix{Float64} = zeros((2 * n + 1)^2, (2 * n + 1)^2)
    for k in 1:4*n^2
        for i in 1:4
            for j in 1:4
                global_index_1::Int64 = get_global_index_node(i, k, n)
                global_index_2::Int64 = get_global_index_node(j, k, n)
                H[global_index_1, global_index_2] += get_v_element(i, j, k, len) + get_t_element(i, j)
                S[global_index_1, global_index_2] += get_s_element(i, j, len)
            end
        end
    end
    return H, S
end

"Handle edge nodes: remove rows and columns, insert diagonal values."
function handle_removing_things(H::Matrix{Float64}, S::Matrix{Float64}, i::Int64)::Nothing
    for j in 1:(2*N+1)^2
        H[i, j] = 0
        H[j, i] = 0
        S[i, j] = 0
        S[j, i] = 0
    end
    S[i, i] = 1
    H[i, i] = -bitwaPodGrunwaldem
end

"Find all global nodes' indexes that are on edge of a box. "
function get_edge_indexes(n::Int16=N)::Int64[]
    output = []
    i = 1
    while i <= 2 * n + 1
        push!(output, i)
        i += 1
    end
    while i <= (2 * n + 1)^2 - (2 * n + 1)
        push!(output, i)
        i += 2 * n
        push!(output, i)
        i += 1
    end
    while i <= (2 * n + 1)^2
        push!(output, i)
        i += 1
    end
    return output
end

# --------------------- ---------------------
function print_table(tab::Matrix{Float64}, n::Int16)::Nothing
    println()
    for i in 1:n
        for j in 1:n
            #  print(tab[i, j], "\t")
            @printf("%.1f\t", tab[i, j])
        end
        println()
    end
    println()
end

function save_matrix_to_file(A::Matrix{Float64}, size::Float16, filename::String)::Nothing
    open("output/$filename", "w") do io
        for i in 1:size
            for j in 1:size
                write(io, "$(A[i, j])\t")
            end
            write(io, "\n")
        end
    end
end

