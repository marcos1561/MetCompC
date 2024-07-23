using Plots
using Random: default_rng, seed!
using Statistics
using DataFrames, CSV

function initialize(N)
    global E, M

    rand_ising2d(m, n=m) = rand(Int8[-1, 1], m, n)

    s0 = rand_ising2d(N);
    s = copy(s0)

    M = mean(s) 
    E = sum( [J * s[i,j] * (s[ifelse(i == N, 1, i+1),j] + 
    s[ifelse(i == 1, N, i-1),j]+ s[i,ifelse(j == N, 1, j+1)] + 
    s[i,ifelse(j == 1, N, j-1)]) for i in 1:N for j in 1:N])

    return s0, s, E, M
end

function animate_gif(s, β, time = 500)
    global E, M
    M_values = []
    M = mean(s)

    anim = @animate for t in 1:time
        rng = default_rng()
        E,M = ising2d_ifelse!(s, β, rng)

        push!(M_values, M)

        p1 = heatmap(s; color=palette([:red,:blue], 3), 
        aspect_ratio=1, size=(800, 800),
        title=string("\nβ = $(β)  Time = $(t)   M = $(round(M,digits=2))"))

        p2 = plot(M_values, title="Magnetization over time", 
        xlabel="Time", ylabel="Magnetization", linewidth=3)

        plot(p1, p2, layout = (2, 1))
    end

    gif(anim, "anim_n$(size(s)[1])_beta$(β).gif", fps = 15)
end
   

function ising2d_ifelse!(s, β, rng=default_rng())
    global E, M
    
    m, n = size(s)
    N2 = m*n
    E_min = -4
    E_max = 4
    prob = [exp(-2*β*E) for E in E_min:E_max]


    for iter in 1:N
        i,j = [rand(1:n) for _ in 1:2]
        
        up = s[ifelse(i == 1, m, i-1), j]
        down = s[ifelse(i == m, 1, i+1), j]
        left = s[i, ifelse(j == 1, n, j-1)]
        right = s[i, ifelse(j == n, 1, j+1)]
        summ = (up + down + left + right)
        prob_idx = summ * s[i,j] - E_min + 1
        dE = 2 * J * summ * s[i,j]

        if rand(rng) < prob[prob_idx]
            s[i,j] = -s[i,j]
            M += 2*s[i,j]/N2
            E += dE
        end
    end

    return E, M
end




const β_crit = log(1+sqrt(2))/2

seed!(4649)

J = 1
β = 0.8
N = 50



### Uncomment to animate
β = 0.2
s,s0,E,M = initialize(N)
animate_gif(s, β)

