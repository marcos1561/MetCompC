@kwdef struct DynCfg 
    g::Float64
    gamma::Float64
    beta::Float64
end

@kwdef struct IntCfg 
    tf::Float64
    dt::Float64
end

function integrate(z0, v0, dyn_cfg, int_cfg)
    dt, tf = int_cfg.dt, int_cfg.tf
    num_steps = trunc(Int, tf/dt) + 1

    z = Array{Float64}(undef, num_steps)
    v = Array{Float64}(undef, num_steps)
    times = Array{Float64}(undef, num_steps)
    
    g, gamma, beta = dyn_cfg.g, dyn_cfg.gamma, dyn_cfg.beta

    z[1] = z0
    v[1] = v0
    times[1] = 0
    for i in 2:num_steps
        dw = randn() * dt^.5
        
        z[i] = z[i-1] + v[i-1] * dt
        
        v1 = v[i-1]
        v[i] = v1 + (g - gamma * v1)*dt + beta * dw
        
        times[i] = times[i-1] + int_cfg.dt
    end

    return times, v, z
end

function stratonovich_mean(z0, v0, dyn_cfg, int_cfg)
    num_sims = 1000
    zs = Array{Float64}(undef, num_sims)
    vs = Array{Float64}(undef, num_sims)

    for i in 1:num_sims
        times, v, z = integrate(z0, v0, dyn_cfg, int_cfg)
        
        zs[i] = z[end]
        vs[i] = v[end]
    end

    # print(zs)
    print(zs[100])
    histogram2d(zs, vs)
end

function v_arg_deriv(j, i, configs, p)
    return (configs.gamma * j - configs.g) * p[j, i]
end

function z_arg_deriv(j, i, configs, p)
    return -j * p[j, i]
end
function v2_arg_deriv(j, i, configs, p)
    return configs.beta^2 * p[j, i]
end

function fokker_planck(dz, dv, num_points, dyn_cfg, int_cfg)
    dt, tf = int_cfg.dt, int_cfg.tf

    p = zeros(Float64, num_points, num_points)
    
    middle = trunc(Int, num_points/2)
    p[middle, middle] = 1
    
    num_steps = trunc(Int, tf/dt) + 1
    times = zeros(Float64, num_steps)
    for i in 2:(num_steps-1)
        old_p = copy(p)
        for i in 2:(num_points-1) # z
            for j in 2:(num_points-1) # v
                df1 = v_arg_deriv(j+1, i, dyn_cfg, old_p) - v_arg_deriv(j-1, i, dyn_cfg, old_p)     
                df1_dv = df1 / 2.0
                
                df2 = z_arg_deriv(j, i+1, dyn_cfg, old_p) - z_arg_deriv(j, i-1, dyn_cfg, old_p)     
                df2_dz = df2 / 2.0
                
                d2f3_dv2 = v2_arg_deriv(j+1, i, dyn_cfg, old_p) + v2_arg_deriv(j-1, i, dyn_cfg, old_p) - 2.0*v2_arg_deriv(j, i, dyn_cfg, old_p)
                
                p[j, i] += (df1_dv + df2_dz + 0.5 * d2f3_dv2) * dt

                # if isnan(x[j])
                #     println("ERRRO")
                # end
            end
        end
        times[i] = times[i-1] + dt
    end

    # length = (num_points - 1)*dx
    
    # a = length / (num_points-1)
    # b = length * (1/2 - num_points/(num_points-1))

    # real_x = range(1, num_points) .* a .+ b

    # return real_x, x ./ dx
    return p
end

using Plots

z0, v0 = 0.0, 0.0

dyn_cfg = DynCfg(
    g = 1,
    gamma = 0.01,
    beta = 0.01,
)

int_cfg = IntCfg(
    tf = 30,
    dt = 0.001,
)


stratonovich_mean(z0, v0, dyn_cfg, int_cfg)

# p = fokker_planck(0.1, 0.1, 200, dyn_cfg, int_cfg)
# cut = 1
# heatmap(p[cut:end-cut, cut:end-cut])

