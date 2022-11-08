include("const.jl")
include("utils.jl")
using Plots
using LinearAlgebra
using Dates

function save_elements()
    open("output/elements.csv", "w") do io
        write(io, "element,local_node,global_node,global_index,local_index,coordinates_x,coordinates_y\n")
        for element in 1:4*N^2
            for local_node in 1:4
                global_node = get_global_index_node(local_node, element)
                global_index = get_global_index_node(local_node, element, N)
                local_index = get_local_index_node(global_index, element, N)
                x, y = get_coordinates_nodes(global_index)
                x *= L0
                y *= L0
                write(io, "$element,$local_node,$global_node,$global_index,$local_index,$x,$y\n")
            end
        end

    end
end

function print_element(k)
    println("element $k")
    for local_node in 1:4
        global_node = get_global_index_node(local_node, k)
        println("local node: $local_node, global node: $global_node")
        print_node(global_node, k)
        println()
    end
end

function print_elements()

    for element_idx in 1:4*N^2
        print_element(element_idx)
        println()
    end
end

function get_psi_nodes(x, y, m=M, ω=OMEGA)
    return exp(-m * ω / 2 * (x^2 + y^2))
end

function get_psi_teo()
    x_size = 20 * 2 * N + 1
    y_size = 20 * 2 * N + 1
    ψ = zeros(x_size, y_size)
    for i in 1:x_size
        for j in 1:y_size
            ψ[i, j] = get_psi_nodes(i * 0.1 * L / (4 * N) - L / 2, j * 0.1 * L / (4 * N) - L / 2)
        end
    end
    return ψ
end

function get_psi(k, ψ, ξ1, ξ2)
    output = 0
    for i in 1:4
        idx = get_global_index_node(i, k)
        row, col = find_row_col_node(idx)
        local_ψ = ψ[(row-1)*20+1, (col-1)*20+1]
        output += local_ψ * g(i, ξ1, ξ2)
    end
    return output
end

function populate_nodes(ψ::Matrix{Float64})
    for i in 1:(2*N+1)^2
        row, col = find_row_col_node(i)
        x, y = get_coordinates_nodes(i)
        ψ[(row-1)*20+1, (col-1)*20+1] = get_psi_nodes(x, y)
    end
end


function calculate_psi_element(element, ψ)
    ξ1 = -1:0.1:1
    ξ2 = -1:0.1:1
    idx = get_global_index_node(1, element)
    row, col = find_row_col_node(idx)
    for i in 1:20
        for j in 1:20
            if i == 1 && j == 1
                continue
            end

            ψ[(row-1)*20+j, (col-1)*20+i] = get_psi(element, ψ, ξ1[i], ξ2[j])
        end
    end
end

function plot_psi(ψ, label)
    x = -L/2:0.1*L/(4*N):L/2
    y = -L/2:0.1*L/(4*N):L/2
    heatmap(x, y, ψ, aspect_ratio=:equal, color=:viridis)
    savefig("output/$label")
end

function ex2()
    ψ = zeros(20 * 2 * N + 1, 20 * 2 * N + 1)
    populate_nodes(ψ)

    for element in 1:4*N^2
        calculate_psi_element(element, ψ)
    end
    plot_psi(ψ, "psi")
    plot_psi(get_psi_teo(), "psi_teo")
end

function solve_eigenproblem(l=L)
    println("current l: $l")
    H, S = get_global_matrix(l)
    E, _ = eigen(H, S)
    return E
end

function main()
    # if output doesnt exist, create it
    if !isdir("output")
        mkdir("output")
    end
    # print_elements()
    # save_elements()

    # ex2()

    # # ex3
    # println(get_matrix(get_s_element))

    # # ex4
    # println(get_matrix(get_t_element))
    H, S = get_global_matrix()
    # print_table(S, (2 * N + 1)^2)
    # save_matrix_to_file(S, (2 * N + 1)^2, "S_initial")
    # save_matrix_to_file(H, (2 * N + 1)^2, "H_initial")


    edges = get_edge_indexes()

    for edge in edges
        handle_removing_things(H, S, edge)
    end

    # print_table(S, (2 * N + 1)^2)
    # save_matrix_to_file(S, (2 * N + 1)^2, "S_final")
    # save_matrix_to_file(H, (2 * N + 1)^2, "H_final")

    # write S to a file

    # E, c = eigen(H, S)
    # println(E)

    matrix = zeros(4)

    for i in 1:4
        matrix[i] = get_v_element(i, i, 7) / get_s_element(i, i)
    end
    plot()
    # println(matrix * R)
    L_num = 20
    L_start = 20
    L_end = 200
    L_step = (L_end - L_start) / L_num
    E = zeros(10, L_num + 1)
    L_arr = L_start:L_step:L_end
    for (index, l) in enumerate(L_arr)
        for j in 1:10
            E[j, index] = solve_eigenproblem(l)[j]
        end
    end
    for i in 1:10
        plot!(L_arr, E[i, :], label="E$i")
        println(E[i, :])
    end
    println(E)

    savefig("output/eigenvalues.png")



end




main()