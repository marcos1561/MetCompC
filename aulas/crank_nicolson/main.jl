using MetCompC, Plots

l = 20
dt = 0.1
k = 0.4
s = 1
tf = 100

f = zeros(l)

source = Vector{Float32}(1:l)
@. source = exp(-(source - l/3)^2 / (10))

plot(f, label="Inicial")
Methods.diffusion_crank!(f, k, dt, s, source, tf)
plot!(f, label="Final")