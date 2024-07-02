using MetCompC, Printf

dynamic_cfg = Molecular.Configs.DynamicCfg(
    ko=2000, ro=0.5, ra=0.5*1.1   
    )
    
int_cfg = Molecular.Configs.IntCfg(dt=0.001)

tf = 4
temp = 1
n_array = [4,9,16,25,36,49,64,81,100]
    
function run(system, tf)
    time = 0
    while time < tf
        Molecular.step!(system)
        time += int_cfg.dt
    end
end

function get_info(n_array, tf)
    times = []
    memory = []
    for n in n_array
        square_n = trunc(Int, n^.5)
        l = dynamic_cfg.ro * square_n
        space_cfg = Molecular.Configs.SpaceCfg(
            length=l, height=l 
        )

        state = Molecular.InitStates.regular_grid(dynamic_cfg.ro/2, temp, space_cfg)
        system = Molecular.System(state, space_cfg, dynamic_cfg, int_cfg)

        info = @timed run(system, tf)
        push!(times, info.time)
        push!(memory, info.bytes * 1e-6)
    end

    return times, memory
end

using Plots

times, memory = get_info(n_array, tf)
scatter(times, label="tempo")
scatter(memory, label="memória")
