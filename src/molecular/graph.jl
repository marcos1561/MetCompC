module Graph

using GLMakie, Printf
using ..Integration: System, State, step!, calc_diff_and_dist!
using ..Quantities

@kwdef struct GraphCfg
    fps::Int
    num_steps_per_frame::Int

    circle_rel::Int = 20
end

function get_circle_base(radius, resolution=20)
    theta = LinRange(0, 2π, resolution+1)

    circle = Vector{Point2f}()
    for id in 1:length(theta)-1
        t1, t2 = theta[id], theta[id+1]
        push!(circle, Point2f(radius * cos(t1), radius * sin(t1))) 
        push!(circle, Point2f(radius * cos(t2), radius * sin(t2))) 
    end

    return circle
end

function update_circles!(circle_points, circle_base, state::State{T}) where {T}
    num_edges = length(circle_base)
    num_p = length(state.x)
    for i in 1:num_p
        begin_id = num_edges * (i-1) + 1
        end_id = begin_id + num_edges - 1 
        center = [state.x[i], state.y[i]]
        circle_points[begin_id:end_id] .= circle_base .+ [center]
    end
end

function animate(system::System{T}, cfg::GraphCfg) where {T}
    x_obs = Observable(system.state.x)
    y_obs = Observable(system.state.y)
    pos_obs = [x_obs, y_obs]
    
    fig, ax = scatter(pos_obs...)
    ax.aspect = 1
    
    # Bordas
    l, h = system.space_cfg.length, system.space_cfg.height
    lines!(ax, [0, l], [0, 0], color=:black)
    lines!(ax, [0, l], [h, h], color=:black)
    lines!(ax, [0, 0], [0, h], color=:black)
    lines!(ax, [l, l], [0, h], color=:black)

    # Círculos
    circle_base = get_circle_base(system.dynamic_cfg.ro/2, cfg.circle_rel)
    circle_points = Vector{Point2f}(undef, system.num_p * length(circle_base))
    
    update_circles!(circle_points, circle_base, system.state)
    circle_points_obs = Observable(circle_points)
    circles_lines = linesegments!(ax, circle_points_obs, color=:black)

    calc_diff_and_dist!(system)
    init_energy = Quantities.energy(system)

    display(fig)
    time = 0.0
    while events(fig).window_open[] 
        for _ in 1:cfg.num_steps_per_frame
            step!(system)
            time += system.int_cfg.dt
        end
        
        energy = Quantities.energy(system)
        temp = Quantities.temperature(system)
        ax.title = @sprintf("E=%.3f | T=%.3f | t=%.3f",energy-init_energy,temp, time)

        update_circles!(circle_points, circle_base, system.state)
        notify(circle_points_obs)

        notify(pos_obs[1])
        notify(pos_obs[2])
        sleep(1/cfg.fps)
    end
end

end

#--
# Animação utilizando o Plots
#--
# @gif for i in 1:num_frames
#     for i in 1:num_steps_per_frame
#         Molecular.step!(system)
#         global time += int_cfg.dt
#     end

#     energy = Molecular.Quantities.energy(system)

#     title=@sprintf("E=%.3f | t=%.3f",energy,time)

#     scatter(particles, system.state.x, system.state.y,
#         title=title, c="red")
# end fps=30