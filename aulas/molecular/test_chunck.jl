include("../../src/MetCompC.jl")

function init_state(num_p_x, num_p_y, offset, radius)
    x = Vector{Float64}()
    y = Vector{Float64}()
    current_y = radius * (offset + 1)
    for i in 1:num_p_y
        current_x = -radius
        for j in 1:num_p_x
            push!(x, current_x + radius * (2 + offset))
            push!(y, current_y)
            current_x = x[end]
        end
        current_y += radius * (2 + offset)
    end
    return x, y
end

num_p_x = 40
num_p_y = 20
num_p = num_p_x * num_p_y
offset = 0.4

dynamic_cfg = MetCompC.Molecular.Configs.DynamicCfg(
    ko=2000, ro=0.5, ra=0.5*1.1   
)

x, y = init_state(num_p_x, num_p_y, offset, dynamic_cfg.ro/2)
space_cfg = MetCompC.Molecular.Configs.SpaceCfg(
    length=maximum(x)+dynamic_cfg.ro/2*(1+offset),
    height=maximum(y)+dynamic_cfg.ro/2*(1+offset),
)

state = MetCompC.Molecular.State{Float64}(
        x, y, 
        (rand(num_p)*2 .- 1), 
        (rand(num_p)*2 .- 1), 
)

chunck = MetCompC.Molecular.Integration.Chuncks(
    trunc(Int, num_p_x/2), trunc(Int, num_p_y/2), space_cfg, state, dynamic_cfg.ro/2,
)

int_cfg = MetCompC.Molecular.Configs.IntCfg(
    dt=0.001
)

graph_cfg = MetCompC.Molecular.Graph.GraphCfg(
    num_steps_per_frame=20,
    fps=60,
    circle_rel=40,
)
    
system = MetCompC.Molecular.System(state, space_cfg, dynamic_cfg, int_cfg, chunck)


MetCompC.Molecular.Graph.animate(system, graph_cfg)