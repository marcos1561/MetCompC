using MetCompC, Printf

dynamic_cfg = Molecular.Configs.DynamicCfg(
    ko=2000, ro=0.5, ra=0.5*1.1   
)

d = dynamic_cfg.ro
space_cfg = Molecular.Configs.SpaceCfg(
    length=d*4,
    height=d*4
    # length=4,
    # height=4
)

int_cfg = Molecular.Configs.IntCfg(
    dt=0.001
)

graph_cfg = Molecular.Graph.GraphCfg(
    num_steps_per_frame=10,
    fps=30,
    circle_rel=40,
)

temp = 1
num_p = 10

tf = 1.

# state = Molecular.InitStates.RandomPos(temp, num_p, space_cfg, dynamic_cfg.ro/2)
state = Molecular.InitStates.regular_grid(dynamic_cfg.ro/2, temp, space_cfg)
# state = Molecular.State{Float64}(
    #     [1., 3], [1, 1], [3, -3], [0, 0]
    # )
    
system = Molecular.System(state, space_cfg, dynamic_cfg, int_cfg)
println(system.num_p)

Molecular.calc_diff_and_dist!(system)
println("Energia K: ", Molecular.Quantities.kinetic_energy(system.state))
println("Energia U: ", Molecular.Quantities.potencial_energy(system))
println("Energia E: ", Molecular.Quantities.energy(system))

Molecular.Graph.animate(system, graph_cfg)

# function run(tf, system)
#     time = 0
#     while time < tf
#         Molecular.step!(system)
#         time += system.int_cfg.dt
#     end
# end
# run(tf, system)
# println(system.num_p)