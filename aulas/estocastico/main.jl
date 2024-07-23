@kwdef struct DynamicCfg 
    alpha::Float64
    beta::Float64
end

@kwdef struct IntCfg 
    dt::Float64
end

function force(x, dynamic_cfg::DynamicCfg)
    alpha = dynamic_cfg.alpha
    return -x^3 + alpha * x
end

function step(x, dynamic_cfg, int_cfg::IntCfg)
    beta = dynamic_cfg.beta 
    x += force(x, dynamic_cfg) * int_cfg.dt + beta * randn() * int_cfg.dt^.5
end

function integrate(tf, xo, dynamic_cfg, int_cfg)
    num_steps = trunc(Int, tf / int_cfg.dt)
    x = Vector{Float64}(undef, num_steps+1)
    x[1] = xo
    for i in 1:num_steps
        x[i+1] = step(x[i], dynamic_cfg, int_cfg)
    end
    return x
end

tf = 1000
xo = 0

int_cfg = IntCfg(
    dt=0.1/2
)

dynamic_cfg = DynamicCfg(
    alpha=0.7,
    beta=0.3,
)

x = integrate(tf, xo, dynamic_cfg, int_cfg)

using Plots

plot(x)

