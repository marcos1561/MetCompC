using Plots

function gauss_seidel(f, source, alpha, tol=0.1, check_freq=4, is_step=false)
    num_lines = size(f)[1]
    num_cols = size(f)[2]
    count = 0
    last_f = copy(f) 
    while true
        for i in 2:num_lines-1
            for j in 2:num_cols-1
                sum_neighboors = f[i-1, j] + f[i+1, j] + f[i, j-1] + f[i, j+1]
                f[i, j] = -alpha * f[i, j] + (1 + alpha)/4 * (sum_neighboors + source[i, j])
            end
        end    

        if is_step == true  
            return 1
        end

        count += 1
        
        if count % check_freq == 0
            max_diff = maximum(abs.(f .- last_f))
            
            if max_diff < tol
                return count
            end

            last_f = copy(f)
        end
    end
end

function apply_bc!(f)
    f[:, 1] .= 1
    f[:, end] .= 0

    l1 = 0.1*2
    l2 = 0.05*2

    @. f[1, :] = exp(-l1*x)
    @. f[end, :] = exp(-l2*x)
end

l = 50
alpha = 0.888*1
num_steps = 1000
tol = 1e-5

x = y = Vector{Float32}(1:l)

f_zeros = zeros(l, l)
f_ones = ones(l, l)
f_rand = rand(l, l)

source = zeros(l, l)

mid_id = trunc(Int, l/2)
source[mid_id, mid_id] = 1

apply_bc!(f_zeros)
apply_bc!(f_ones)
apply_bc!(f_rand)

# f[:, 1] .= 1
# f[:, end] .= 1
# f[1, :] .= 1
# f[end, :] .= 1

num_steps_zeros = gauss_seidel(f_zeros, source, alpha, tol)
num_steps_ones = gauss_seidel(f_ones, source, alpha, tol)
num_steps_rand = gauss_seidel(f_rand, source, alpha, tol)

p1 = heatmap(f_zeros, title="zeros | $(num_steps_zeros)")
p2 = heatmap(f_ones, title="ones | $(num_steps_ones)")
p3 = heatmap(f_rand, title="rand | $(num_steps_rand)")
# plot(p1, p2, p3, layout=(3,1))
# display(p1)
savefig(p1, "por_do_sol_deitado.png")


# num_frames = 60
# heatmap(f_zeros)
# @gif for i in 1:num_frames
#     gauss_seidel(f_zeros, source, alpha, is_step=true)
#     heatmap!(f_zeros)
# end
