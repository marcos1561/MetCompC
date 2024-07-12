include("../../src/MetCompC.jl")

space_cfg = MetCompC.Molecular.Configs.SpaceCfg(
    length=10,
    height=20,
    # length=4,
    # height=4
)

num_p = 5
state = MetCompC.Molecular.State{Float64}(
        (rand(num_p) .+ 1) * 4, 
        (rand(num_p) .+ 1) * 4, 
        (rand(num_p) .+ 1) * 4, 
        (rand(num_p) .+ 1) * 4, 
    )

chunck = MetCompC.Molecular.Integration.Chuncks(
    5, 5, space_cfg, state, 0.5,
)

dynamic_cfg = MetCompC.Molecular.Configs.DynamicCfg(
    ko=2000, ro=0.5, ra=0.5*1.1   
)

int_cfg = MetCompC.Molecular.Configs.IntCfg(
    dt=0.001
)

graph_cfg = MetCompC.Molecular.Graph.GraphCfg(
    num_steps_per_frame=10,
    fps=30,
    circle_rel=40,
)
    
system = MetCompC.Molecular.System(state, space_cfg, dynamic_cfg, int_cfg, chunck)


MetCompC.Molecular.Graph.animate(system, graph_cfg)