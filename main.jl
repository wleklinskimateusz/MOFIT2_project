include("const.jl")
include("utils.jl")
include("matrix.jl")
include("net.jl")
include("psi.jl")
using Plots
using LinearAlgebra
using Dates

function save_elements(n::Int16, l::Float64)::Nothing
    open("output/elements.csv", "w") do io
        write(io, "element,local_node,global_node,global_index,local_index,coordinates_x,coordinates_y\n")
        for element in 1:4*n^2
            for local_node in 1:4
                global_node = get_global_index_node(local_node, element, n)
                global_index = get_global_index_node(local_node, element, n)
                local_index = get_local_index_node(global_index, element, n)
                x, y = get_coordinates_nodes(global_index, l, n)
                x *= L0
                y *= L0
                write(io, "$element,$local_node,$global_node,$global_index,$local_index,$x,$y\n")
            end
        end

    end
end


function ex2(n::Int16)::Matrix{Float64}
    l::Float64 = 100 / L0
    ψ::Matrix{Float64} = zeros(20 * 2 * n + 1, 20 * 2 * n + 1)
    populate_nodes(ψ, get_psi_nodes, n, l)

    for element in 1:4*n^2
        calculate_psi_element(element, ψ, n)
    end
    plot_psi(ψ, "psi", n, l)
    plot_psi(get_psi_teo(l, n), "psi_teo", n, l)
    return ψ
end

function ex4a()::Nothing
    n::Int16 = 2
    l::Float64 = 100 / L0
    elements = [6, 7, 10, 11]
    println()
    for element in elements
        println("Element: $element")
        for i in 1:4
            x, y = get_coordinates_nodes(get_global_index_node(i, element, n), l, n)
            value = get_v_matrix_element(i, i, 7, l, n) / get_s_element(i, i, l, n) * R
            println("$(x*L0), $(y*L0):\t $value")
        end
        println("")
    end
end

function ex5(n::Int16)::Nothing
    states = 7
    plot()
    L_num = 20
    L_start = 20
    L_end = 200
    L_step = (L_end - L_start) / L_num
    E = zeros(states, L_num + 1)
    L_arr::Vector{Float64} = L_start:L_step:L_end
    L_arr ./= L0
    for (index::Int64, l::Float64) in enumerate(L_arr)
        E[:, index], _ = solve_eigenproblem(l, n, states)
    end
    for i in 1:states # loop through eigenstates
        plot!(L_arr * L0, E[i, :] * R, label="E$i")
    end
    xlabel!("L [nm]")
    ylabel!("E [eV]")
    savefig("output/eigenvalues_$n.png")
    return nothing
end

function ex6(l::Float64, n::Int16)::Tuple{Vector{Float64},Matrix{Float64}}
    E, c = solve_eigenproblem(l, n, 6)
    ψ = zeros(20 * 2 * n + 1, 20 * 2 * n + 1)
    for state in 1:6
        populate_nodes(ψ, (i, n, l) -> c[state, i], n, l)
        for element in 1:4*n^2
            calculate_psi_element(element, ψ, n)
        end
        plot_psi(ψ, "psi_$state", n, l)
    end
    return E, c
end

function get_normalised_c(c::Vector{Float64}, S::Matrix{Float64})::Vector{Float64}
    return c ./ sqrt(c' * S * c)
end

function get_x_analytical(E1::Float64, E2::Float64, tmax::Int64, A::Float64=1.0)::Vector{Float64}
    T = 2 * π / (E2 - E1)
    println(T)
    t = 1:Δt:(tmax*Δt)
    return A * cos.(T * t)
end

function ex7(E::Vector{Float64}, c::Matrix{Float64}, l::Float64, n::Int16, tmax::Int64=100)::Nothing
    H, S = get_global_matrix(l, n)
    X = get_global_x_matrix(l, n)

    c1::Vector{Float64} = get_normalised_c(c[1, :], S)
    c2::Vector{Float64} = get_normalised_c(c[2, :], S)
    d::Vector{ComplexF64} = c1 + c2

    E1, E2 = E[1:2]
    x_theo = get_x_analytical(E1, E2, tmax, 0.6)
    x::Vector{Float64} = zeros(tmax)
    ψ = zeros(20 * 2 * n + 1, 20 * 2 * n + 1)
    anim::Animation = Animation()

    times = 1:Δt:(tmax*Δt)
    for t in 1:tmax
        d = inv(get_matrix_left(H, S, n)) * get_matrix_right(H, S, n) * d
        # hermitian transpose
        x[t] = real(d' * X * d)
        d_abs2 = real(d' .* d)
        populate_nodes(ψ, (i, n, l) -> d_abs2[i], n, l)
        for element in 1:4*n^2
            calculate_psi_element(element, ψ, n)
        end
        h = plot_psi(ψ, "t=$t", n, l, false)
        frame(anim, h)
    end
    gif(anim, "output/psi.gif")

    plot(times, x_theo, label="x_theo")
    plot!(times, x, label="x")
    xlabel!("t")
    ylabel!("x")
    savefig("output/x.png")
    return nothing
end



function main()
    # if output doesnt exist, create it
    if !isdir("output")
        mkdir("output")
    end
    n::Int16 = 2
    l::Float64 = 80 / L0
    # @time print_elements(n, l)
    # save_elements(n, l)

    # @time ex2(n)

    # ex3
    # size::Int16 = 4
    # println("Lokalna macierz S.")
    # print_table(get_matrix((i, j) -> get_s_element(i, j, l, n)), size)
    # print_table(get_matrix((i, j) -> get_s_element(i, j, l, n)) / get_a(l, n)^2 * 36, size)

    # ex4
    # println("Lokalna macierz t.")
    # print_table(get_matrix(get_t_element), size)
    # print_table(get_matrix(get_t_element) * 6 * 2 * M, size)

    # @time ex4a()

    # @time ex5(n)

    # n = 10
    E, c = @time ex6(100.0, n)
    # solve_eigenproblem(l, n, 6)
    # @time ex7(E, c, 100.0, n)

end




main()