function delta_x(x, y, wo, beta, dw, dt)
    (-wo * y - beta^2 / 2 * x) * dt - beta * y * dw
end

function delta_y(x, y, wo, beta, dw, dt)
    (wo * x - beta^2 / 2 * y) * dt + beta * x * dw
end

function integrate(;tf, dt, wo, beta)
    xs = [1.]
    ys = [0.]
    
    xi = [1.]
    yi = [0.]
    
    xe = [1.]
    ye = [0.]
    
    times = [0.]

    num_steps = trunc(Int, tf/dt)
    w = 0
    for j in 1:num_steps
        dw = dt^.5 * rand()
        w += dw
        fase = wo * dt * j + beta * w

        
        # Strat
        xj, yj = xs[end], ys[end]
        c = (wo * dt + beta * dw) / 2
        d = 1 + c^2
        x = (xj - c*yj)/d
        y = (yj + c*xj)/d
        push!(xs, 2*x - xj)
        push!(ys, 2*y - yj)

        # Ito
        xj, yj = xi[end], yi[end]
        push!(xi, xj + delta_x(xj, yj, wo, beta, dw, dt))
        push!(yi, yj + delta_y(xj, yj, wo, beta, dw, dt))
        
        # Exata
        push!(xe, cos(fase))
        push!(ye, sin(fase))
        
        push!(times, times[end] + dt)
    end

    return times, xi, yi, xs, ys, xe ,ye
end

using Plots

times, xi, yi, xs, ys, xe, ye = integrate(
    tf = 1,
    dt = 0.01,
    wo = 1,
    beta = 1.5,
)

p1 = plot(times, [xi xs xe], label=["ito" "strat" "exata"], ylabel="x")
p2 = plot(times, [yi ys ye], label=["ito" "strat" "exata"], ylabel="y")
plot(p1, p2, layout=(2, 1), xlabel="tempo")