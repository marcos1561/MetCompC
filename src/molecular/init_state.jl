module InitStates

using Distributions
using ..Integration: State

function RandomPos(temp, num_particles, space_cfg, r)
    maxwell_d = Normal(0, sqrt(temp))

    state = State{Float64}(
        x = rand(num_particles) * (space_cfg.length - 2*r) .+ r,
        y = rand(num_particles) * (space_cfg.height - 2*r) .+ r,
        vx = rand(maxwell_d, num_particles),
        vy = rand(maxwell_d, num_particles),
    )

    return state
end

function regular_grid(radius, temp, space_cfg)
    y = radius
    x_array = Vector{Float64}()
    y_array = Vector{Float64}()
    while y < space_cfg.height - radius
        x = radius
        while x < space_cfg.length - radius
            push!(x_array, x)
            push!(y_array, y)
            x += 2*radius
        end
        y += 2*radius
    end
    
    maxwell_d = Normal(0, sqrt(temp))
    num_particles = length(x_array)
    state = State(
        x=x_array, y=y_array,
        vx = rand(maxwell_d, num_particles),
        vy = rand(maxwell_d, num_particles),
    )

    return state
end

end