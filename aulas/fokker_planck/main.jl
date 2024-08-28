@kwdef struct DynCfg 
    k::Float64
    gamma::Float64
    beta::Float64
end

@kwdef struct IntCfg 
    tf::Float64
    dt::Float64
end

function integrate(x0, dyn_cfg, int_cfg)
    dt, tf = int_cfg.dt, int_cfg.tf
    num_steps = trunc(Int, tf/dt) + 1

    x = Array{Float64}(undef, num_steps)
    times = Array{Float64}(undef, num_steps)
    
    k, gamma, beta = dyn_cfg.k, dyn_cfg.gamma, dyn_cfg.beta

    k_gamma = k / gamma

    dt2 = int_cfg.dt^2

    x[1] = x0
    times[1] = 0
    for i in 2:num_steps
        dw = randn() * dt2
        xi = x[i-1]
        x[i] = xi - k_gamma * xi * int_cfg.dt + beta * dw
        times[i] = times[i-1] + int_cfg.dt
    end

    return times, x
end

function fokker_planck(dx, num_points, dyn_cfg, int_cfg)
    dt, tf = int_cfg.dt, int_cfg.tf
    k, gamma, beta = dyn_cfg.k, dyn_cfg.gamma, dyn_cfg.beta

    x = zeros(Float64, num_points)
    
    middle = trunc(Int, num_points/2)
    x[middle] = 1
    
    num_steps = trunc(Int, tf/dt) + 1
    times = zeros(Float64, num_steps)
    for i in 2:(num_steps-1)
        old_x = copy(x)
        for j in 2:(num_points-1)
            f_j_plus = (convert(Float64, j)+1. - convert(Float64, middle)) * old_x[j+1]
            f_j_minus = (convert(Float64, j)-1. - convert(Float64, middle)) * old_x[j-1]

            deriv = (f_j_plus - f_j_minus)/2.0
            deriv2 = old_x[j+1] - 2*old_x[j] + old_x[j-1]
            
            dp_dt = k/gamma * deriv + 1/2 * (beta/gamma)^2  * deriv2
            x[j] = old_x[j] + dp_dt * dt

            if isnan(x[j])
                println("ERRRO")
            end
        end
        times[i] = times[i-1] + dt
    end

    length = (num_points - 1)*dx
    
    a = length / (num_points-1)
    b = length * (1/2 - num_points/(num_points-1))

    real_x = range(1, num_points) .* a .+ b

    return real_x, x ./ dx
end

using Plots

dyn_cfg = DynCfg(
    k = 1,
    gamma = 1,
    beta = 0.1,
)

int_cfg = IntCfg(
    tf = 3.8,
    dt = 0.01,
)

x0 = 0

# times, x = integrate(x0, dyn_cfg, int_cfg)

x, p = fokker_planck(0.01, 500, dyn_cfg, int_cfg)

println(p[50])
# println(Vector(x))

plot(x, p)
# plot(x)


