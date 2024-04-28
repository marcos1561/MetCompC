function sem_ponto(f)
    for i in 1:100
        f[2:end-1] = f[1:end-2] + f[3:end] 
    end
end

function com_ponto(f)
    for i in 1:100
        @. f[2:end-1] = f[1:end-2] + f[3:end] 
    end
end

l = 10000000
f1 = Vector{Float32}(1:l)
f2 = Vector{Float32}(1:l)

@time sem_ponto(f1)
@time sem_ponto(f1)
println("======")
@time com_ponto(f2)
@time com_ponto(f2)

difff = zeros(l)
@. difff = (f2 - f1)^2

# println(f1)
# println(f2)
println(sum(difff))
