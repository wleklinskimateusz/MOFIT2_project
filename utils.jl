function get_length_atomic(l::Float64)::Float64
    return l / L0
end

function get_length_metric(l::Float64)::Float64
    return l * L0
end

function get_a(l::Float64, n::Int16)::Float64
    return l / (2 * n)
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



