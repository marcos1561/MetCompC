using Distributions

@kwdef struct Pos{T}
    x::Vector{T}
    y::Vector{T}
    vx::Vector{T}
    vy::Vector{T}
end

@kwdef struct SpaceCfg
    length::Float32
    height::Float32
end

@kwdef struct ForcePars 
    n::Int
    b::Float32

    ko::Float32
    ro::Float32
    ra::Float32
end

function calc_temp(pos)
    n = size(pos.vx)[1]
    return sum(pos.vx.^2 + pos.vy.^2)/(2*n)
end

function harmonic_pot(r, f_pars::ForcePars)
    if r > f_pars.ra
        return 0.0
    end

    return -2*f_pars.ko * (0.5*f_pars.ra^2 - f_pars.ra*f_pars.ro - 0.5*r^2 + f_pars.ro * r)
end

function calc_potencial(pos, diff, f_pars::ForcePars)
    ue = 0
    n = size(pos.x)[1]

    for i in 1:n
        for j in (i+1):n
            dist = (diff[1, i, j]^2 + diff[2, i, j]^2)^.5  
            ue += harmonic_pot(dist, f_pars)
        end
    end 
    return ue        
end

function calc_kinetic(pos)
    return sum((pos.vx.^2 + pos.vy.^2) ./ 2.0)
end

function calc_energy(pos, diff, f_pars)
    # return calc_potencial(pos, diff, f_pars)
    return calc_kinetic(pos) + calc_potencial(pos, diff, f_pars)
end

function calc_force(pos, f_pars::ForcePars)
    n = size(pos.x)[1]

    diff = Array{Float32, 3}(undef, 2, n, n)
    for j in 1:n
        for i in (j+1):n
            diff[1, i, j] = pos.x[i] - pos.x[j]
            diff[2, i, j] = pos.y[i] - pos.y[j]
            diff[1, j, i] = pos.x[j] - pos.x[i]
            diff[2, j, i] = pos.y[j] - pos.y[i]
        end
    end

    println(calc_energy(pos, diff, f_pars))
    
    nf = f_pars.n
    b = f_pars.b
    forces = zeros(Float32, 2, n)
    for j in 1:n
        for i in (j+1):n
            dist = (diff[1, i, j]^2 + diff[2, i, j]^2)^.5        
            if dist == 0
                println("Erro: r nulo!")
            end

            # fo = 2 * nf / dist^(2*nf + 2) - b*nf/dist^(nf)
            
            fo = 0
            if dist < f_pars.ra
                fo = -f_pars.ko * (dist - f_pars.ro)
            end

            fx = fo * diff[1, i, j] / dist
            fy = fo * diff[2, i, j] / dist

            forces[1, i] +=  fx
            forces[2, i] +=  fy
            
            forces[1, j] += -fx
            forces[2, j] += -fy
        end
    end

    return forces
end

function step!(pos, f_pars, space, dt)
    n = size(pos.x)[1]

    forces1 = calc_force(pos, f_pars)
    for i in 1:n
        pos.x[i] += dt * pos.vx[i]
        pos.y[i] += dt * pos.vy[i]
    end
    
    forces2 = calc_force(pos, f_pars)
    forces_mean = 1/2 .* (forces1 .+ forces2)
    for i in 1:n
        pos.vx[i] += dt * forces_mean[1, i]
        pos.vy[i] += dt * forces_mean[2, i]
    end

    for i in 1:n
        if (pos.x[i] > space.length) || (pos.x[i] < 0)
            pos.vx[i] *= -1
        end
        if (pos.y[i] > space.height) || (pos.y[i] < 0)
            pos.vy[i] *= -1
        end
    end
end




n = 25

dt = 0.001
tf = 1

space_cfg = SpaceCfg(
    length = 2.2,
    height = 2.2,
)
f_pars = ForcePars(
    n=4, b=1,
    ko=200, ro=0.5, ra=0.5*1.1,
)

temp = 1


maxwell_d = Normal(0, sqrt(temp))
pos = Pos{Float32}(
    x =rand(n) * space_cfg.length,
    y = rand(n) * space_cfg.height,
    vx = rand(maxwell_d, n),
    vy = rand(maxwell_d, n),
)

using Plots
@gif for i in 1:400
    for i in 1:2
        step!(pos, f_pars, space_cfg, dt)
    end

    scatter(pos.x, pos.y, legend=false)
    xlims!(0, space_cfg.length)
    ylims!(0, space_cfg.height)
end fps=10


