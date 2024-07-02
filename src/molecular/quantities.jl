module Quantities

export energy, kinetic_energy, potencial_energy, harmonic_potencial, temperature

function kinetic_energy(state)
    return 0.5 * sum(state.vx.^2 + state.vy.^2)
end

function harmonic_potencial(r, pars)
    if r > pars.ra
        return 0.0
    end
    
    b = pars.ra - pars.ro
    return pars.ko/2. * ((r - pars.ro)^2 - b^2) 
end

function potencial_energy(system)
    ue = 0

    for i in 1:system.num_p
        for j in (i+1):system.num_p
            dist = system.dists[i, j]
            ue += harmonic_potencial(dist, system.dynamic_cfg)
        end
    end 
    return ue 
end

function energy(system)
    return kinetic_energy(system.state) + potencial_energy(system)
end
   
function temperature(system)
    return kinetic_energy(system.state) / system.num_p
end

end
