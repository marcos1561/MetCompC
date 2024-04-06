using MetCompC, Plots

l = 30
k = 0.1
tf = 4
dt = 0.01

sigma = 20
f = Vector{Float32}(1:l)
@. f = exp(-(f - l/2)^2/sigma)

plot(f, label="Inicial")
Methods.deriva_ftcs!(f, k, tf, dt)
plot!(f, label="Final")
