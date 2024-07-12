module Integration

export State, System, calc_diff_and_dist!, calc_forces!, step!

using ..Configs

@kwdef struct State{T}
    x::Vector{T}
    y::Vector{T}
    vx::Vector{T}
    vy::Vector{T}
end

struct System{T}
    state::State{T}
    space_cfg::SpaceCfg
    dynamic_cfg::DynamicCfg
    int_cfg::IntCfg

    diffs::Array{Float64, 3}
    forces::Array{Float64, 2}
    dists::Array{Float64, 2}

    num_p::Int

    function System{T}(state, space_cfg, dynamic_cfg, int_cfg) where T
        num_p = length(state.x)

        diffs = Array{Float32, 3}(undef, 2, num_p, num_p)
        forces = Array{Float32, 2}(undef, 2, num_p)
        dists = zeros(Float32, num_p, num_p)
        
        new(state, space_cfg, dynamic_cfg, int_cfg, diffs, forces, dists, num_p)
    end
end
System(state::State{T}, space_cfg, dynamic_cfg, int_cfg) where {T} = System{T}(state, space_cfg, dynamic_cfg, int_cfg)

function calc_diff_and_dist!(system::System{T}) where {T}
    state = system.state
    diffs = system.diffs
    dists = system.dists

    for j in 1:system.num_p
        for i in (j+1):system.num_p
            dx = state.x[i] - state.x[j]
            dy = state.y[i] - state.y[j]

            diffs[1, i, j] = dx
            diffs[2, i, j] = dy
            diffs[1, j, i] = -dx
            diffs[2, j, i] = -dy
            
            r = (dx^2 + dy^2)^.5 
            dists[i, j] = dists[j, i] = r
        end
    end
end

function calc_forces!(system::System{T}) where {T}
    diffs = system.diffs
    dists = system.dists
    forces = system.forces
    f_pars = system.dynamic_cfg
    n = system.num_p

    forces .= 0.0

    calc_diff_and_dist!(system)
    # n = length(state.x)
    # for j in 1:n
    #     for i in (j+1):n
    #         diffs[1, i, j] = state.x[i] - state.x[j]
    #         diffs[2, i, j] = state.y[i] - state.y[j]
    #         diffs[1, j, i] = state.x[j] - state.x[i]
    #         diffs[2, j, i] = state.y[j] - state.y[i]
    #     end
    # end

    # info.energy = calc_energy(pos, diffs, f_pars)
    
    for j in 1:n
        for i in (j+1):n
            # dist = (diffs[1, i, j]^2 + diffs[2, i, j]^2)^.5        
            # system.dists[i, j] = dist
            # system.dists[j, i] = dist
            
            dist = dists[i, j]
            if dist == 0
                println("Erro: r nulo!")
            end

            # fo = 2 * nf / dist^(2*nf + 2) - b*nf/dist^(nf)
            
            fo = 0
            if dist < f_pars.ra
                fo = -f_pars.ko * (dist - f_pars.ro)
            end

            fx = fo * diffs[1, i, j] / dist
            fy = fo * diffs[2, i, j] / dist

            forces[1, i] +=  fx
            forces[2, i] +=  fy
            
            forces[1, j] -= fx
            forces[2, j] -= fy
        end
    end
end

function step!(system::System{T}) where {T}
    state = system.state
    dt = system.int_cfg.dt

    n = length(system.state.x)

    info = @timed calc_forces!(system)
    # println("Force time: ", info.time)
    # println("Force mem: ", info.bytes)
    for i in 1:n
        state.x[i] += dt * state.vx[i] + 0.5*dt^2 * system.forces[1, i] 
        state.y[i] += dt * state.vy[i] + 0.5*dt^2 * system.forces[2, i]
    end
    
    forces_last = copy(system.forces)
    calc_forces!(system)
    forces_mean = 1/2 .* (forces_last .+ system.forces)
    for i in 1:n
        state.vx[i] += dt * forces_mean[1, i]
        state.vy[i] += dt * forces_mean[2, i]
    end

    space_cfg = system.space_cfg
    r = system.dynamic_cfg.ro/2
    for i in 1:n
        if ((state.x[i]+r) > space_cfg.length) || ((state.x[i]-r) < 0)
            state.vx[i] *= -1.
        end
        if ((state.y[i]+r) > space_cfg.height) || ((state.y[i]-r) < 0)
            state.vy[i] *= -1.
        end
    end
end

end