function delta_x(x, y, wo, beta, dw, dt)
    (-wo * y - beta^2 / 2 * x) * dt - beta * y * dw
end

function delta_y(x, y, wo, beta, dw, dt)
    (wo * x - beta^2 / 2 * y) * dt + beta * x * dw
end

function integrate(;tf, dt, wo, beta)
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

        xj, yj = xi[end], yi[end]

        push!(xi, xj + delta_x(xj, yj, wo, beta, dw, dt))
        push!(yi, yj + delta_y(xj, yj, wo, beta, dw, dt))
        push!(xe, cos(fase))
        push!(ye, sin(fase))
        push!(times, times[end] + dt)
    end

    return times, xi, yi, xe ,ye
end

using Plots

times, xi, yi, xe, ye = integrate(
    tf = 3,
    dt = 0.01,
    wo = 1,
    beta = 0.5,
)

p1 = plot(times, [xi xe], label=["ito" "exata"])
p2 = plot(times, [yi ye], label=["ito" "exata"])
plot(p1, p2, layout=(2 , 1))