using Plots
using Random,Distributions
using LinearAlgebra
using Printf

function inicia()
    sqN=10
    N=sqN^2
    n=2
    dt=0.001
    tmax=2.
    T=1.0
    dy = 0.1
    
    r0=0.5
    ddx = r0
    ddy = r0
    X = ddx*sqN
    Y = ddy*sqN
    
    ra = 1.1*r0
    k0 = 2000.0
    
    c=distinguishable_colors(N)
    x = Array{Float64,1}(undef,N)
    y = Array{Float64,1}(undef,N)
    for i = 1:N
        j=(i-1)%sqN;l=div((i-1),sqN)
        x[i]=j*ddx+ddx/2;y[i]=l*ddy+ddy/2
    end
    d=Normal(0,sqrt(T))
    Random.seed!(1234)
    vx=rand(d,N)
    vy=rand(d,N)
    Tm=sum(vx .^2+vy .^2)/(2N)
    println("Temperatura=$Tm")
    
    # N=2
    # c=distinguishable_colors(N)
    # x = [1., 1.01]
    # y = [2., 2.]
    # vx = [3., -3.]
    # vy = [0.0, 0.0]
    # X, Y = 5., 5.

    return x,y,vx,vy,Tm,N,X,Y,dy,c,dt,tmax,n,r0,k0,ra
end

function plot_frame(x,y,Tm,t,N,X,Y,dy,c,V,K,TE)
	 title=@sprintf("T=%.2f t=%.2f TE=%.3f V=%.2f K=%.2f",Tm,t,TE,V,K)
	 global p=scatter(x,y,title=title,c=c,ms=7,legend=false)
	[annotate!(x, y+dy, Plots.text(string(i), 12)) for (i,x,y) in zip(1:N,x,y)]
	 plot!(p,xlimits=(0,X),ylimits=(0,Y))
	 return p
end


function diff(x,y)
	 diff_x = x .- x'
	 diff_y = y .- y'
	 diff_r = sqrt.(diff_x .* diff_x + diff_y .* diff_y)
	 diff_r += Diagonal(ones(N))
	 return diff_r,diff_x,diff_y
end

function contorno(x,y,vx,vy,X,Y)
	 zx=findall(t->(t<0)||(t>X),x)
	 vx[zx]=-vx[zx]
	 zy=findall(t->(t<0)||(t>Y),y)
	 vy[zy]=-vy[zy]
return vx,vy
end

function force_mod(r)
    #f = 12 .*(1 ./r .^(2*n+2) -1 ./r .^(n+2))
    f = -k0 .*(r .- r0)
return f
end

function force(x,y)
    diff_r,diff_x,diff_y=diff(x,y)
    f_mod=force_mod.(diff_r)
    f_mod[diagind(f_mod)] .= 0.0
    z=findall(t->(t > ra),diff_r)
    f_mod[z] .= 0.0
    fxx = f_mod .* diff_x ./ diff_r
    fyy = f_mod .* diff_y ./ diff_r
    fx = sum(fxx,dims=2)
    fy = sum(fyy,dims=2)
    return fx,fy
end

function potential(r)
    Vr =  0.5 .*k0 .*(r .^2 .- ra .^2) .- k0 .*r0 .*(r .- ra)
return Vr
end

function Total_Energy(x,y,vx,vy)
    V = PotentialEnergy(x,y)
    K = Kinetic(vx,vy)
    TE = V + K
    return V,K,TE
end
function Kinetic(vx,vy)
        K = sum(vx .^2+vy .^2)/2
    return K
end

function PotentialEnergy(x,y)
    diff_r,diff_x,diff_y=diff(x,y)
    
    # println(size(diff_r))
    # println(diff_r)
    
    
    Vmatrix = potential(diff_r)
    
    
    Vmatrix[diagind(Vmatrix)] .= 0.0
    z=findall(t->(t > ra),diff_r)
    Vmatrix[z] .= 0.0
    
    # println(size(Vmatrix))
    # println(Vmatrix)
    
    # V=0.25*sum(Vmatrix)
    V=0.5*sum(Vmatrix)
    
    return V
end

function move(x,y,vx,vy,fxi,fyi,dt)
    mdt2=dt^2/2; h=dt/2
    x += vx*dt + mdt2*fxi 
    y += vy*dt + mdt2*fyi
    fxf,fyf=force(x,y)
    vx += h*(fxi+fxf)
    vy += h*(fyi+fyf) 
    return x,y,vx,vy,fxf,fyf
end

function evol(tmax,x,y,vx,vy,fx,fy,dt)
    # t=0;icount=0
    # num_steps_per_frame = 10

    # num_steps = trunc(Int, tmax/dt)

    # num_frames = num_steps / num_steps_per_frame
    # @gif for _ in 1:num_frames 
    #     for _ in 1:num_steps_per_frame
    #         t +=dt;icount += 1
    #         x,y,vx,vy,fx,fy=move(x,y,vx,vy,fx,fy,dt)
    #         vx,vy=contorno(x,y,vx,vy,X,Y)
    #     end
    #     V,K,TE = Total_Energy(x,y,vx,vy)

    #     plot_frame(x,y,Tm, fx[1] ,N,X,Y,dy,c,V,K,TE)
    # end

    t=0;icount=0
    while t < tmax
        t +=dt;icount += 1
        x,y,vx,vy,fx,fy=move(x,y,vx,vy,fx,fy,dt)
        vx,vy=contorno(x,y,vx,vy,X,Y)
        # if mod(icount,10) == 0
        #     V,K,TE = Total_Energy(x,y,vx,vy)
        #     p=plot_frame(x,y,Tm,t,N,X,Y,dy,c,V,K,TE)
        #     display(p)
        #     sleep(.1)
        # end
    end
end
x,y,vx,vy,Tm,N,X,Y,dy,c,dt,tmax,n,r0,k0,ra=inicia()
fx,fy=force(x,y)

# V,K,TE = Total_Energy(x,y,vx,vy)
# println("V=$(V),K=$(K),E=$(TE)")

println(N)

@time evol(tmax,x,y,vx,vy,fx,fy,dt)
#readline()