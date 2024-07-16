module Integration

export State, System, calc_diff_and_dist!, calc_forces!, step!

using ..Configs

@kwdef struct State{T}
    x::Vector{T}
    y::Vector{T}
    vx::Vector{T}
    vy::Vector{T}
end

struct Chuncks{T}
    num_cols::Int
    num_rows::Int
    chunck_length::Float64
    chunck_height::Float64
    space_cfg::SpaceCfg
    steps_to_update::Int

    state::State{T}
    
    neighbors::Matrix{Vector{CartesianIndex{2}}}
    
    chunck_particles::Array{Int, 3}
    num_particles_in_chunck::Array{Int, 2}
end
function Chuncks(num_cols, num_rows, space_cfg, state, particle_r)
    neighbors = Matrix{Vector{CartesianIndex{2}}}(undef, num_rows, num_cols)
    for i in 1:(num_rows-1)
        neighbors[i, 1] = [
            CartesianIndex(i+1, 1),
            CartesianIndex(i+1, 2),
            CartesianIndex(i, 2),
        ]
        
        for j in 2:(num_cols-1)
            neighbors[i, j] = [
                CartesianIndex(i+1, j-1),
                CartesianIndex(i+1, j),
                CartesianIndex(i+1, j+1),
                CartesianIndex(i, j+1),
            ]
        end

        neighbors[i, num_cols] = [
            CartesianIndex(i+1, num_cols-1),
            CartesianIndex(i+1, num_cols),
        ]
    end
    
    for j in 1:(num_cols-1)
        neighbors[num_rows, j] = [CartesianIndex(num_rows, j+1)]
    end
    neighbors[num_rows, num_cols] = []
    
    chunck_length = space_cfg.length / num_cols
    chunck_height = space_cfg.height / num_rows

    nc = (ceil(0.5*chunck_length/particle_r) + 1) * ((ceil(0.5*chunck_height/particle_r) + 1))
    nc = trunc(Int, ceil(nc*1.1))
    chunck_particles = Array{Int}(undef, nc, num_rows, num_cols)
    num_particles_in_chunck = zeros(Int, num_rows, num_cols)

    Chuncks(num_cols, num_rows, chunck_length, chunck_height, space_cfg, 100, 
        state, neighbors, chunck_particles, num_particles_in_chunck)
end

mutable struct Info
    time::Float64
    time_step::Int
end

struct System{T}
    state::State{T}
    space_cfg::SpaceCfg
    dynamic_cfg::DynamicCfg
    int_cfg::IntCfg

    chuncks::Chuncks

    diffs::Array{Float64, 3}
    forces::Array{Float64, 2}
    dists::Array{Float64, 2}

    num_p::Int

    info::Info
    function System{T}(state, space_cfg, dynamic_cfg, int_cfg, chuncks) where T
        num_p = length(state.x)

        diffs = Array{Float32, 3}(undef, 2, num_p, num_p)
        forces = Array{Float32, 2}(undef, 2, num_p)
        dists = zeros(Float32, num_p, num_p)
        
        new(state, space_cfg, dynamic_cfg, int_cfg, chuncks, diffs, forces, dists, num_p, Info(0, 0))
    end
end
System(state::State{T}, space_cfg, dynamic_cfg, int_cfg, chuncks) where {T} = System{T}(
    state, space_cfg, dynamic_cfg, int_cfg, chuncks)
    
function update_chuncks!(chuncks)
    state = chuncks.state

    space_h = chuncks.space_cfg.height
    chunck_l, chunck_h = chuncks.chunck_length, chuncks.chunck_height

    chuncks.num_particles_in_chunck .= 0

    num_p = length(state.x)
    for i in 1:num_p
        row_id = trunc(Int, div(-state.y[i] + space_h, chunck_h)) + 1
        col_id = trunc(Int, div(state.x[i], chunck_l)) + 1
        
        # println(state.x[i], "|", state.y[i])
        # println(row_id, "| ", col_id)
        # println("========")

        row_id -= row_id == (chuncks.num_rows + 1) ? 1 : 0
        col_id -= col_id == (chuncks.num_cols + 1) ? 1 : 0

        p_i = chuncks.num_particles_in_chunck[row_id, col_id] + 1
        chuncks.chunck_particles[p_i, row_id, col_id] = i 
        chuncks.num_particles_in_chunck[row_id, col_id] += 1
    end
end

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

function calc_interaction!(system, i, j)
    state = system.state
    forces = system.forces
    f_pars = system.dynamic_cfg

    dx = state.x[i] - state.x[j]
    dy = state.y[i] - state.y[j]
    dist = (dx^2 + dy^2)^.5 

    if dist == 0
        println("Erro: r nulo!")
    end

    if dist > f_pars.ra
        return
    end

    fo = -f_pars.ko * (dist - f_pars.ro)
    fx = fo * dx / dist
    fy = fo * dy / dist

    forces[1, i] +=  fx
    forces[2, i] +=  fy
    
    forces[1, j] -= fx
    forces[2, j] -= fy
end

function calc_forces!(system::System{T}) where {T}
    system.forces .= 0
    
    chuncks = system.chuncks
    for col in 1:chuncks.num_cols
        for row in 1:chuncks.num_rows
            np = chuncks.num_particles_in_chunck[row, col]
            chunck = @view chuncks.chunck_particles[1:np, row, col]
            neighbors = chuncks.neighbors[row, col]
            
            for i in 1:np
                p1_id = chunck[i]
                for j in (i+1):np
                    # println("$(i), $(j), $(np)")
                    calc_interaction!(system, p1_id, chunck[j])
                end
                
                for neighbor_id in neighbors
                    nei_np = chuncks.num_particles_in_chunck[neighbor_id]
                    nei_chunck = @view chuncks.chunck_particles[1:nei_np, neighbor_id]
                    
                    for j in 1:nei_np
                        calc_interaction!(system, p1_id, nei_chunck[j])
                    end
                end
            end
        end
    end
end

function calc_forces_normal!(system::System) 
    system.forces .= 0
    for i in 1:system.num_p
        for j in (i+1):system.num_p
            calc_interaction!(system, i, j)
        end
    end
end

function step!(system::System{T}) where {T}
    state = system.state
    dt = system.int_cfg.dt

    n = length(system.state.x)

    if system.info.time_step % system.chuncks.steps_to_update == 0
        update_chuncks!(system.chuncks)
    end

    # calc_diff_and_dist!(system)
    # calc_forces_normal!(system)
    calc_forces!(system)
    
    # info = @timed calc_forces!(system)
    # println("Force time: ", info.time)
    # println("Force mem: ", info.bytes)
    for i in 1:n
        state.x[i] += dt * state.vx[i] + 0.5*dt^2 * system.forces[1, i] 
        state.y[i] += dt * state.vy[i] + 0.5*dt^2 * system.forces[2, i]
    end
    
    forces_last = copy(system.forces)
    # calc_diff_and_dist!(system)
    # calc_forces_normal!(system)
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

    system.info.time += system.int_cfg.dt
    system.info.time_step += 1
end

end