using LinearAlgebra
using Plots

function get_psi_teo(l::Float64, n::Int16)::Matrix{Float64}
    x_size = 20 * 2 * n + 1
    y_size = 20 * 2 * n + 1
    ψ = zeros(x_size, y_size)
    for i in 1:x_size
        for j in 1:y_size
            ψ[i, j] = get_psi_nodes(i * 0.1 * l / (4 * n) - l / 2, j * 0.1 * l / (4 * n) - l / 2)
        end
    end
    return ψ
end

"Compute wave function values only in nodes."
function get_psi_nodes(i::Int64, n::Int16, l::Float64, m::Float64=M, ω::Float64=OMEGA)::Float64
    x, y = get_coordinates_nodes(i, l, n)
    return exp(-m * ω / 2 * (x^2 + y^2))
end

"Compute wave function values only in nodes."
function get_psi_nodes(x::Float64, y::Float64, m::Float64=M, ω::Float64=OMEGA)::Float64
    return exp(-m * ω / 2 * (x^2 + y^2))
end

"Compute wave function values in between nodes in coordinates given in frame of reference by ξ1, ξ2."
function get_psi(k::Int64, ψ::Matrix{Float64}, ξ1::Float64, ξ2::Float64, n::Int16)::Float64
    output::Float64 = 0
    for i in 1:4
        idx::Int64 = get_global_index_node(i, k, n)
        row, col = find_row_col_node(idx, n)
        local_ψ = ψ[(row-1)*20+1, (col-1)*20+1]
        output += local_ψ * g(i, ξ1, ξ2)
    end
    return output
end

"Fill values in ψ matrix in places corresponding to nodes."
function populate_nodes(ψ::Matrix{Float64}, populate_node::Function, n::Int16, l::Float64)
    for i in 1:(2*n+1)^2
        row, col = find_row_col_node(i, n)
        ψ[(row-1)*20+1, (col-1)*20+1] = populate_node(i, n, l)
    end
end


function calculate_psi_element(element::Int64, ψ::Matrix{Float64}, n::Int16)::Nothing
    ξ1 = -1:0.1:1
    ξ2 = -1:0.1:1
    idx = get_global_index_node(1, element, n)
    row, col = find_row_col_node(idx, n)
    for i in 1:20
        for j in 1:20
            if i == 1 && j == 1
                continue
            end

            ψ[(row-1)*20+j, (col-1)*20+i] = get_psi(element, ψ, ξ1[i], ξ2[j], n)
        end
    end
end

function plot_psi(ψ::Matrix{Float64}, label::String, n::Int16, l::Float64, save::Bool=true)::Plots.Plot{Plots.GRBackend}
    x = -l/2:0.1*l/(4*n):l/2
    y = -l/2:0.1*l/(4*n):l/2
    h = heatmap(x, y, ψ, aspect_ratio=:equal, color=:viridis)
    if save
        savefig("output/$label")
    end
    return h
end

function solve_eigenproblem(len::Float64, n::Int16, states::Int64)::Tuple{Vector{Float64},Matrix{Float64}}
    H, S = get_global_matrix(len, n)
    E, c = eigen(H, S)
    i = 1
    while E[i] < 0
        i += 1
    end
    return E[i:i+states-1], c[i:i+states-1, :]
end

