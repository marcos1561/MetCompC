square_func(w, dt, dw) = w^2 * dw + w * dt

function integrate(t_max, dt_small, dt_big, integrand)
    num_steps_big = trunc(Int, t_max/dt_big)
    
    w_small = 0
    w_big = 0
    t = 0
    
    times = []
    w_list = []
    int_small = []
    int_big = []
    int_big_int = 0
    
    int_ep = 0
    for j in 1:num_steps_big
        dw_big = 0
        for i in 1:10
            dw_small = dt_small^.5 * randn()
            dw_big += dw_small
            int_ep += integrand(w_small, dt_small, dw_small)
            w_small += dw_small 
        end
        push!(int_small, int_ep)
        int_big_int += integrand(w_big, dt_big, dw_big)
        push!(int_big, int_big_int)
        
        w_big += dw_big
        t += dt_big

        push!(w_list, w_big)     
        push!(times, t)
    end

    int_ito = 1/3 * w_list.^3
    return times, int_small, int_big, int_ito
end

using Plots

t_max = 100
dt_small = 0.1
dt_big = 10 * dt_small

times, small, big, ito = integrate(t_max, dt_small, dt_big, square_func)

plot(times, [small, big, ito], label=["Pequeno" "Grande" "Ito"])

