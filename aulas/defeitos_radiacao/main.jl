using Plots

l = 30
d = 3
s = 0.5

x = Vector{Float32}(1:l)

f = zeros(l)

sigma = 1
source = zeros(l)
@. source = exp(-(x-l/2)^2 / (2*sigma^2))
source[2] += d * f[1]
source[l-1] += d * f[l-1]

M = zeros(l-2, l-2)

M[1, 1] = -2*d -s
M[1, 2] = d
M[end, end-1] = d
M[end, end] = -2*d -s
for i in 2:l-3
    M[i, i] = -2*d - s
    M[i, i+1] = d
    M[i, i-1] = d
end

M_inv = inv(M)

plot(f, label="inicial")
f[2:end-1] = -M_inv * source[2:end-1]
plot!(f, label="final")
plot!(source, label="fonte")

# gui()

