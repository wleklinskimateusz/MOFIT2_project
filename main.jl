include("utils.jl")
include("const.jl")
using Plots

function print_elements()
    for element_idx in 1:4*N^2
        println("element $element_idx")
        for local_node in 1:4
            global_node = get_global_index_node(local_node, element_idx)
            println("local node: $local_node, global node: $global_node")
            print_node(global_node, element_idx)
            println()
        end
    end
end

function get_psi_nodes(x, y, m=M, omega=OMEGA)
    return exp(-m * omega / 2 * (x^2 + y^2))
end

function get_psi_teo()
    x_size = 20 * 2 * N + 1
    y_size = 20 * 2 * N + 1
    psi = zeros(x_size, y_size)
    for i in 1:x_size
        for j in 1:y_size
            psi[i, j] = get_psi_nodes(i * 0.1 * L / (4 * N) - L / 2, j * 0.1 * L / (4 * N) - L / 2)
        end
    end
    return psi
end

function get_psi(k, psi, ksi_1, ksi_2)
    output = 0
    for i in 1:4
        idx = get_global_index_node(i, k)
        row, col = find_row_col_node(idx)
        local_psi = psi[(row-1)*20+1, (col-1)*20+1]
        output += local_psi * g(i, ksi_1, ksi_2)
    end
    return output
end

function populate_nodes(psi::Matrix{Float64})
    for i in 1:(2*N+1)^2
        row, col = find_row_col_node(i)
        x, y = get_coordinates_nodes(i)
        psi[(row-1)*20+1, (col-1)*20+1] = get_psi_nodes(x, y)
    end
end


function calculate_psi_element(element, psi)
    ksi_1 = -1:0.1:1
    ksi_2 = -1:0.1:1
    idx = get_global_index_node(1, element)
    row, col = find_row_col_node(idx)
    for i in 1:20
        for j in 1:20
            if i == 1 && j == 1
                continue
            end

            psi[(row-1)*20+j, (col-1)*20+i] = get_psi(element, psi, ksi_1[i], ksi_2[j])
        end
    end
end

function plot_psi(psi, label)
    x = -L/2:0.1*L/(4*N):L/2
    y = -L/2:0.1*L/(4*N):L/2
    heatmap(x, y, psi, aspect_ratio=:equal, color=:viridis)
    savefig("output/$label")

end

function main()
    # if output doesnt exist, create it
    if !isdir("output")
        mkdir("output")
    end
    # print_elements()

    psi = zeros(20 * 2 * N + 1, 20 * 2 * N + 1)
    populate_nodes(psi)

    for element in 1:4*N^2
        calculate_psi_element(element, psi)
    end
    plot_psi(psi, "psi")
    plot_psi(get_psi_teo(), "psi_teo")
end


main()