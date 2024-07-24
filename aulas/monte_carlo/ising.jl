@kwdef struct Pars
    j::Float64
    temperature::Float64
    kb::Float64
end

@kwdef struct System{T}
    spins::Matrix{T}
    pars::Pars
end

@kwdef struct PreCalcVars{T}
    directions::Vector{CartesianIndex{2}}
    probs_dict::Dict{Int8, Float64}
    linear_to_cartesian_idx::T
end
PreCalcVars(system) = PreCalcVars(
    directions = [CartesianIndex(0, 1), CartesianIndex(1, 0), CartesianIndex(0, -1), CartesianIndex(-1, 0)],
    probs_dict = possibles_delta_energys(system.pars),
    linear_to_cartesian_idx = CartesianIndices(system.spins),
)

function possibles_delta_energys(pars::Pars)
    values = Vector{Tuple{Int8, Float64}}()
    for i in 1:9
        de = -5 + i
        if de <= 0
            push!(values, (-5 + i, 1.0))
        else
            push!(values, (-5 + i, exp(-2.0*pars.j*de / (pars.kb * pars.temperature))))
        end
    end

    return Dict(values)
end

 
function sum_neighbors_energy(center_id, spins, directions)
    num_rows, num_cols = size(spins)
    sum = 0
    for dir in directions
        index = center_id + dir

        if index[1] > num_rows
            index = CartesianIndex(1, index[2])
        end
        if index[1] < 1
            index = CartesianIndex(num_rows, index[2])
        end
        if index[2] > num_cols
            index = CartesianIndex(index[1], 1)
        end
        if index[2] < 1
            index = CartesianIndex(index[1], num_cols)
        end

        sum += spins[index] 
    end

    sum *= spins[center_id]
    return sum
end

function create_random_spins(num_rows, num_cols)
    spins = Matrix{Int8}(undef, num_rows, num_cols)
    indeces = rand(num_rows, num_cols) .> 0.5
    spins[indeces] .= 1
    spins[.!indeces] .= -1
    return spins
end

function step!(system, pre_calc_vars)
    center_id = pre_calc_vars.linear_to_cartesian_idx[rand(1:length(system.spins))]
    energy_k = sum_neighbors_energy(center_id, system.spins, pre_calc_vars.directions)

    accept_flip = false
    if energy_k < 0
        accept_flip = true
    elseif rand() < pre_calc_vars.probs_dict[energy_k]
        accept_flip = true
    end

    if accept_flip
        system.spins[center_id] *= -1
    end
end

function integrate!(system, num_steps, pre_calc_vars=nothing)
    if pre_calc_vars === nothing
        pre_calc_vars = PreCalcVars(system)
    end
    
    for _ in 1:num_steps
        step!(system, pre_calc_vars)
    end
end

function magnetization(spins)
    return sum(spins)/length(spins)
end

function magnetization_curve(;temp_range, pars, size, num_steps)
    mags = []
    for t in temp_range
        current_pars = Pars(
            j = pars.j,
            temperature = t,
            kb = pars.kb,
        )

        system = System(
            spins = create_random_spins(size, size),
            pars = current_pars
        )

        integrate!(system, num_steps)
        push!(mags, magnetization(system.spins))
    end

    return mags
end

temp_range = Vector(0.1:0.1:4)
num_steps = 500000
mags= magnetization_curve(
    temp_range = temp_range, 
    pars = Pars(1, 1, 1),
    size = 40,
    num_steps = num_steps,
)

using Plots, Colors

plot(temp_range, abs.(mags), 
    marker=:circle,
    xlabel="Temperatura",
    ylabel="Magnetização",
    title="Número de passos por simulação: $(num_steps)"
)

# system = System(
#     spins = create_random_spins(40, 40),
#     pars = Pars(
#         j = 1,
#         temperature = 0.1,
#         kb = 1,
#     ),
# )



# cmap = cgrad([:blue, :red], [-1, 1])
# h1 = heatmap(system.spins, color=cmap)
# integrate!(system, 50000)
# println(magnetization(system.spins))
# h2 = heatmap(system.spins, color=cmap)

# plot(h1, h2, layout=(2, 1), size=(800/2, 1200/2))
