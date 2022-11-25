using Printf

function get_s_element(i::Int64, j::Int64, len::Float64, n::Int16)::Float64
    output::Float64 = 0
    for l in 1:3
        for s in 1:3
            output += w[l] * w[s] * g(j, p[l], p[s]) * g(i, p[l], p[s])
        end
    end
    return output * get_a(len, n)^2 / 4
end

function dev_g_ξ_1(j::Int64, x1::Float64, x2::Float64)::Float64
    return (g(j, x1, x2 + ε) - g(j, x1, x2 - ε)) / (2 * ε)
end

function dev_g_ξ_2(j::Int64, x1::Float64, x2::Float64)::Float64
    return (g(j, x1 + ε, x2) - g(j, x1 - ε, x2)) / (2 * ε)
end

function get_t_matrix_element(i::Int64, j::Int64)::Float64
    output = 0
    for l in 1:3
        for n in 1:3
            output += w[l] * w[n] * (dev_g_ξ_1(j, p[l], p[n]) * dev_g_ξ_1(i, p[l], p[n]) + dev_g_ξ_2(j, p[l], p[n]) * dev_g_ξ_2(i, p[l], p[n]))
        end
    end
    return output / (2 * M)
end


function get_v_matrix_element(i::Int64, j::Int64, k::Int64, len::Float64, n::Int16)::Float64
    output = 0
    for l in 1:3
        for s in 1:3
            output += w[l] * w[s] * g(j, p[l], p[s]) * g(i, p[l], p[s]) * (get_x_real(p[l], k, len, n)^2 + get_y_real(p[s], k, len, n)^2)
        end
    end
    return get_a(len, n)^2 / 4 * M * OMEGA^2 / 2 * output
end

# 
function get_x_matrix_element(i::Int64, j::Int64, k::Int64, len::Float64, n::Int16)::Float64
    output = 0
    for l in 1:3
        for s in 1:3
            output += w[l] * w[s] * g(j, p[l], p[s]) * g(i, p[l], p[s]) * get_x_real(p[l], k, len, n)
        end
    end
    return output
end

function get_global_x_matrix(len::Float64, n::Int16)::Matrix{Float64}
    output::Matrix{Float64} = zeros((2 * n + 1)^2, (2 * n + 1)^2)
    for k in 1:4*n^2
        for i in 1:4
            for j in 1:4
                global_index_1::Int64 = get_global_index_node(i, k, n)
                global_index_2::Int64 = get_global_index_node(j, k, n)
                output[global_index_1, global_index_2] += get_x_matrix_element(i, j, k, len, n)
            end
        end
    end
    return output
end

"Create global S and H matrixes."
function get_global_matrix(len::Float64, n::Int16)::Tuple{Matrix{Float64},Matrix{Float64}}
    H::Matrix{Float64} = zeros((2 * n + 1)^2, (2 * n + 1)^2)
    S::Matrix{Float64} = zeros((2 * n + 1)^2, (2 * n + 1)^2)
    for k in 1:4*n^2
        for i in 1:4
            for j in 1:4
                global_index_1::Int64 = get_global_index_node(i, k, n)
                global_index_2::Int64 = get_global_index_node(j, k, n)
                H[global_index_1, global_index_2] += get_v_matrix_element(i, j, k, len, n) + get_t_matrix_element(i, j)
                S[global_index_1, global_index_2] += get_s_element(i, j, len, n)
            end
        end
    end
    for edge in get_edge_indexes(n)
        handle_removing_things(H, S, edge)
    end
    return H, S
end

"Handle edge nodes: remove rows and columns, insert diagonal values."
function handle_removing_things(H::Matrix{Float64}, S::Matrix{Float64}, i::Int64)::Nothing
    H[i, :] .= 0
    H[:, i] .= 0
    S[i, :] .= 0
    S[:, i] .= 0

    S[i, i] = 1
    H[i, i] = -bitwaPodGrunwaldem
    return nothing
end

"Find all global nodes' indexes that are on edge of a box. "
function get_edge_indexes(n::Int16)::Vector{Int64}
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

function get_matrix_left(H::Matrix{Float64}, S::Matrix{Float64}, n::Int16)::Matrix{ComplexF64}
    output::Matrix{ComplexF64} = zeros((2 * n + 1)^2, (2 * n + 1)^2)
    for i in 1:(2*n+1)^2
        for j in 1:(2*n+1)^2
            output[i, j] = S[i, j] - Δt / (2 * im) * H[i, j]
        end
    end
    return output
end

function get_matrix_right(H::Matrix{Float64}, S::Matrix{Float64}, n::Int16)::Matrix{ComplexF64}
    output::Matrix{ComplexF64} = zeros((2 * n + 1)^2, (2 * n + 1)^2)
    for i in 1:(2*n+1)^2
        for j in 1:(2*n+1)^2
            output[i, j] = S[i, j] + Δt / (2 * im) * H[i, j]
        end
    end
    return output
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