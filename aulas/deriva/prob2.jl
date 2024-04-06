using MetCompC, Plots

l = 30
k = 1
tf = 40
dt = 0.1

sigma = 20
f = Vector{Float32}(1:l)
@. f = exp(-(f - l/2)^2/sigma)

p =plot(f, label="Inicial")
Methods.deriva_lax!(f, k, tf, dt)
plot!(f, label="Final")
display(p)

function animate(f, k, tf, dt)
    num_steps = Int(tf/dt)
    t = 0
    @gif for i=1:num_steps
        t_round = round(t, digits=2)
        plot(f, title="t=$(t_round)")
        Methods.deriva_lax_step!(f, k)
        t += dt
    end
end
# animate(f, k, tf, dt)

