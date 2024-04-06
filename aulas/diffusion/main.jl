using MetCompC, Plots

# Configurações
l = 100
tf = 140
dt = 0.1

ks = [0.1, 0.2, 0.3]

# Matrix contendo todos os estados
# A j-ésima coluna é o estado referente ao k[j]
fs = zeros(l, length(ks))
fs[1, :] .= 1

# Gerando a animação
labels = ["k = $(k)" for k in ks]
labels = reshape(labels, (1, length(ks)))

t = 0
num_steps = Int(tf/dt)
@gif for i=1:num_steps
    t_round = round(t, digits=2) 
    plot(fs, label=labels, title="t = $(t_round)")
    
    # Integração ocorre aqui
    for i in eachindex(ks)
        Methods.diffusion_ftcs_step!(@view(fs[:, i]), ks[i])
    end

    global t += dt
end every 10