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

    state::State{T}
    
    # neighbors::Vector{Matrix{Int}}
    neighbors::Array{Int, 4}
    
    chunk_particles::Array{Int, 3}
    num_particles_in_chunck::Array{Int, 2}
end
function Chuncks(num_cols, num_rows, space_cfg, state, particle_r)
    neighbors = Array{Int}(undef, num_rows, num_cols, 2, 4)
    neighbors .= -1
    for i in 1:(num_rows-1)
        neighbors[i, 1, :, 1] = [i+1, 1]
        neighbors[i, 1, :, 2] = [i+1, 2]
        neighbors[i, 1, :, 3] = [i, 2]
        
        for j in 2:(num_cols-1)
            neighbors[i, j, :, 1] = [i+1, j-1]
            neighbors[i, j, :, 2] = [i+1, j]
            neighbors[i, j, :, 3] = [i+1, j+1]
            neighbors[i, j, :, 4] = [i, j+1]
        end
    end

    chunck_length = space_cfg.length / num_cols
    chunck_height = space_cfg.height / num_rows

    area = chunck_height * chunck_length
    nc = trunc(Int, ceil(area / (π * particle_r^2)))
    chunck_particles = Array{Int}(undef, num_rows, num_cols, trunc(Int, ceil(4*nc)))
    num_particles_in_chunck = zeros(Int, num_rows, num_cols)

    Chuncks(num_cols, num_rows, chunck_length, chunck_height, space_cfg, state, neighbors, chunck_particles, num_particles_in_chunck)
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

    function System{T}(state, space_cfg, dynamic_cfg, int_cfg, chuncks) where T
        num_p = length(state.x)

        diffs = Array{Float32, 3}(undef, 2, num_p, num_p)
        forces = Array{Float32, 2}(undef, 2, num_p)
        dists = zeros(Float32, num_p, num_p)
        
        new(state, space_cfg, dynamic_cfg, int_cfg, chuncks, diffs, forces, dists, num_p)
    end
end
System(state::State{T}, space_cfg, dynamic_cfg, int_cfg, chuncks) where {T} = System{T}(state, space_cfg, dynamic_cfg, int_cfg, chuncks)
    
function update_chunck!(chuncks)
    state = chuncks.state

    l, h = chuncks.space_cfg.length, chuncks.space_cfg.height

    chuncks.num_particles_in_chunck .= 0

    num_p = length(state.x)
    for i in 1:num_p
        row_id = trunc(Int, div(-state.y[i] + chuncks.space_cfg.height, chuncks.chunck_height)) + 1
        col_id = trunc(Int, div(state.x[i], chuncks.chunck_length)) + 1
        
        # println(state.x[i], "|", state.y[i])
        # println(row_id, "| ", col_id)
        # println("========")

        row_id -= row_id == chuncks.num_rows ? 1 : 0
        col_id -= col_id == chuncks.num_cols ? 1 : 0

        p_i = chuncks.num_particles_in_chunck[row_id, col_id] + 1
        chuncks.chunk_particles[row_id, col_id, p_i] = i 
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
    diffs = system.diffs
    forces = system.forces
    f_pars = system.dynamic_cfg

    dist = system.dists[i, j]
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

function calc_forces!(system::System{T}) where {T}
    system.forces .= 0
    
    chuncks = system.chuncks
    for col in 1:chuncks.num_cols
        for row in 1:chuncks.num_rows
            np = chuncks.num_particles_in_chunck[row, col]
            chunck = @view chuncks.chunk_particles[row, col, 1:np]
            neighbors = @view chuncks.neighbors[row, col, :, :]
            
            for i in 1:np
                p1_id = chunck[i]
                for j in (i+1):np
                    # println("$(i), $(j), $(np)")
                    calc_interaction!(system, p1_id, chunck[j])
                end
                
                for nei_id in 1:size(neighbors)[2]
                    nei_row = neighbors[1, nei_id]
                    if nei_row == -1
                        continue
                    end
                    nei_col = neighbors[2, nei_id]
                    
                    nei_np = chuncks.num_particles_in_chunck[nei_row, nei_col]
                    nei_chunck = @view chuncks.chunk_particles[nei_row, nei_col, 1:nei_np]
                    
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

    update_chunck!(system.chuncks)
    
    calc_diff_and_dist!(system)
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
    calc_diff_and_dist!(system)
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
end

end