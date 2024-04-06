module Methods

"""
Passo temporal da equação da difusão

    ∂f/∂t = D * ∂²∂/∂x²

utilizando ftcs com bordas fixas.
    
# Argumentos
- `x`: Estado inicial
- `k`: `D*dt/dx²`
"""
function diffusion_ftcs_step!(x, k)
    x[2:end-1] += k * (x[3:end] + x[1:end-2] -2*x[2:end-1])
end

"""
Integra a equação de difusão até o tempo `tf` dado
o estado inicial `x` em `t=0` com o método ftcs. 

Ver `diffusion_ftcs_step!`
para mais detalhes. 
"""
function diffusion_ftcs!(x, tf, dt, k)
    t = 0
    while t < tf
        diffusion_ftcs_step!(x, k)
        t += dt
    end
end

"""
Integra a equação da difusão com fonte

    ∂f/∂t = D * ∂²f/∂x² - S * f + s(x)

com o método Crank-Nicolson com condições de borda periódicas.  

# Argumentos
- f: Estado inicial 
- k: `D * dt / dx²` 
- dt: Tamanho do passo temporal 
- s: Constante da destruição.
- source: fonte em cada x. 
- tf: Tempo final. 
"""
function diffusion_crank!(f, k, dt, s, source, tf)
    l = length(f)
    C = zeros(l, l)
    B = zeros(l, l)
    
    C[1, l] = -k 
    C[1, 1] = (2 + 2*k + s *dt)
    C[1, 2] = -k 
    C[l, l] = (2 + 2*k + s *dt)
    C[l, l-1] = k 
    C[l, 1] = -k 
    
    B[1, l] = k 
    B[1, 1] = (2 - 2*k - s *dt)
    B[1, 2] = k 
    B[l, l] = (2 - 2*k - s *dt)
    B[l, l-1] = k 
    B[l, 1] = k 
    for i in 2:(l-1)
        C[i, i-1] = -k 
        C[i, i] = (2 + 2*k + s *dt)
        C[i, i+1] = -k
        
        B[i, i-1] = k 
        B[i, i] = (2 - 2*k - s *dt)
        B[i, i+1] = k
    end

    C_inv = inv(C)
    
    t = 0
    while t < tf
        f .= C_inv * (B*f + source*dt)
        t += dt
    end
end

function deriva_ftcs_step!(f, k)
    f_copy = copy(f) 
    
    f[1] -= k*(f_copy[2] - f_copy[end])
    f[end] -=  k*(f_copy[1] - f_copy[end-1])

    l = length(f)
    for i in 2:l-1
        f[i] += -k*(f_copy[i+1] - f_copy[i-1])
    end
end

"""
Integra a equação da deriva utilizando o método
ftcs

- OBS: Para essa equção esse método nunca é estável.
"""
function deriva_ftcs!(f, k, tf, dt)
    t = 0
    while t < tf 
        deriva_ftcs_step!(f, k)
        t += dt
    end
end

"""
Passo temporal da equação da deriva

    ∂f/∂t = -v ∂f/∂x

utilizando o método lax, dado a condição inicial `f`
e `k=v*dt/dx`.
    
# Argumentos
- `x`: Estado inicial
- `k`: `D*dt/dx²`
"""
function deriva_lax_step!(f, k)
    f_copy = copy(f) 

    f[1] = 1/2 * (f_copy[end] + f_copy[2]) - k/2 * (f_copy[2] - f_copy[end])
    f[end] = 1/2 * (f_copy[end-1] + f_copy[1]) - k/2 * (f_copy[1] - f_copy[end-1])

    l = length(f)
    for i in 2:l-1
        f[i] = 1/2 * (f_copy[i-1] + f_copy[i+1]) - k/2 * (f_copy[i+1] - f_copy[i-1])
    end
end

"""
Integra a equação da deriva até o tempo `tf`
dado a condição incial `f` em `t=0`.

Ver `deriva_lax_step!` para mais detalhes. 
"""
function deriva_lax!(f, k, tf, dt)
    t = 0
    while t < tf 
        deriva_lax_step!(f, k)
        t += dt
    end
end

end