struct System{T}
    spins::Matrix{T}
end

struct Pars
    j::Float64
    temperature::Float64
    kb::Float64
end

function possible_delta_E(pars::Pars)
    values = []
    for i in 1:8
        de = -5 + i
        if de <= 0
            push!(values, (-5 + i, 1.0))
        else
            push!(values, (-5 + i, exp(-pars.j*de / (pars.kb * pars.temperature))))
        end
    end

    return Dict(values)
end

directions = [CartesianIndex(0, 1), CartesianIndex(1, 0), CartesianIndex(0, -1), CartesianIndex(-1, 0)]
function sum_neighbors(center::CartesianIndex, system::System)
    num_rows, num_cols = size(system.spins)
    sum = 0
    for dir in directions
        index = center + dir

        if index[1] > num_rows
            index = CartesianIndex(index[1], 1)
        end
        if index[1] < 0
            index = CartesianIndex(index[1], num_cols)
        end
        if index[2] > num_cols
            index = CartesianIndex(1, index[2])
        end
        if index[2] < 0
            index = CartesianIndex(num_rows, index[2])
        end
        sum += system.spins[center] 
    end

    sum *= system.spins[center]
    return sum
end

function create_random(num_rows, num_cols)
    spins = Matrix{Int8}(undef, num_rows, num_cols)
    indeces = rand(n, n) .> 0.5
    spins[indeces] .= 1
    spins[.!indeces] .= -1
    return System(spins)
end

function step!(system, pars)
    index = CartesianIndex(rand)

    sum = sum_neighbors()
end

create_random(4)

# pars = Pars(1, 1, 1)
# de_to_prob = possible_delta_E(pars)


# de_to_prob[3]

# function calc_diff_energy(i, j, system::System)
    
# end