using Random 
using GLMakie # for plots
using DelimitedFiles # for importing data
using Printf # for printing out numbers and strings
using StatsBase # for histograms
using LsqFit # for curve fitting



function initialize(N,L) 

    max_rows = Int8(ceil(N/L)) + 1 # placed at a distance of 1

    x = repeat(1:L-1,max_rows)
    x = x[1:N]

    y = repeat(1:max_rows,1,L-1)
    y = y'[:]
    y = y[1:N]
    

    vx = rand(N)
    vx = vx .-sum(vx)

    vy = rand(N)
    vy = vy .-sum(vy)

    return x,y,vx,vy

end

function thermalization(x,y,vx,vy,n_therm,g,L,tau)

    fx,fy = force(x,y,g)
    
    vx +=  tau/2*fx # only go half step in first calculation
    vy +=  tau/2*fy  

    for n in 1:n_therm

        x,y = update_xy(x,y,vx,vy,tau)
        x,y,vx,vy = reflection(x,y,vx,vy,L)
        vx,vy = update_vxvy(x,y,vx,vy,g,tau)

    end


    return x,y,vx,vy
end

function simulation_loop(x,y,vx,vy,n_loop,Y,VX,VY,E,g,L,n_skip,tau,animation)

    counter = 0

    for n in 1:n_loop

        

        x,y = update_xy(x,y,vx,vy,tau)
        x,y,vx,vy = reflection(x,y,vx,vy,L)

        if mod(n,n_skip) == 0 # measure and animate every n_skip steps

            counter += 1
        

            Y[counter,:] = y
            VX[counter,:] = vx
            VY[counter,:] = vy
            E[counter] = energy(x,y,vx,vy,g)

            if animation == true

                animate(x,y,n,L,n_skip)
    
            end

        end

        vx,vy = update_vxvy(x,y,vx,vy,g,tau)

    end

    return Y,VX,VY,E

end

function update_xy(x,y,vx,vy,tau)


    x += tau*vx
    y += tau*vy

    return x,y

end

function reflection(x,y,vx,vy,L)

    ref_left = x .< 0
    ref_right = x .> L
    ref_bottom = y .< 0

    x[ref_left] = -x[ref_left]
    x[ref_right] = -x[ref_right] .+ 2*L
    y[ref_bottom] = -y[ref_bottom]

    vx[ref_left] = -vx[ref_left]
    vx[ref_right] = -vx[ref_right]
    vy[ref_bottom] = -vy[ref_bottom]

    return x,y,vx,vy

end

function update_vxvy(x,y,vx,vy,g,tau)

    fx,fy = force(x,y,g)

    vx += tau*fx
    vy += tau*fy

    return vx,vy

end

function animate(x,y,n,L,n_skip)


    if n == n_skip
        fig = Figure(); display(fig)
        global ax = Axis(fig[1,1])
        points = Point2f[] # initialize so it can be deleted also in the first iteration
        global l = scatter!(ax, points, color = :purple)
    end

    N = length(x)

    points = Observable(Point2f[(x[1],y[1])]) # Point2f[(x[1],y[1])] so that it is of type 1-element Vector{Point{2, Float32}}, not Point           
    
    delete!(ax.scene, l)
    l = scatter!(ax, points, color = :purple)
    limits!(ax, 0, L, 0, L)


    for i = 2:N
        new_point = Point2f(x[i],y[i])
        points[] = push!(points[], new_point)
     
    end
   
    sleep(1/10000)  # refreshes the display!

end

function velocity_histogram(data) 

    h = fit(Histogram, data;)

    binedges = collect(h.edges[1])
    binwidth = abs(binedges[2] - binedges[1])
    numbins = length(binedges) - 1
    bincenters = binedges[1:end-1] .+ 0.5*binwidth


    binheights = 1/binwidth*(h.weights .+ 1)/(length(data) + numbins + 1)
    sigma = (1/binwidth)*sqrt.( (1/(length(data) + numbins + 2))*(binwidth*binheights.*(1 .- binwidth*binheights)) )


    #fig = Figure(); display(fig)
    Axis(fig[1, 1], title="Distribution of velocities")
    

    barplot!(bincenters, binheights, width=binwidth,gap=0,label="MD simulation")
    errorbars!(bincenters, binheights, Vector{Float64}(vec(2*sigma)), color = :red)

    # define model for fit

    m(v,c) = 1/c[1]*v.*exp.(-v.^2 /c[1]/2) # fit parameter(s) have to be of type Array
    c0 = [100.0] # initial guess is important !!!



    fit_func = LsqFit.curve_fit(m,bincenters,binheights,c0)


    x_fit = LinRange(minimum(binedges), maximum(binedges),100)
    y_fit = m(x_fit,fit_func.param)

    plot!(x_fit,y_fit,label="LS-Fit")
    axislegend()


    @printf("kT/m from velocity histogram: %.3f \n", fit_func.param[1])
    @printf("kT/m from the mean velocity: %.3f \n", mean(abs.(data).^2)/2)


    return nothing

end

function density_histogram(data,g) 

    h = fit(Histogram, data;)

    binedges = collect(h.edges[1])
    binwidth = abs(binedges[2] - binedges[1])
    numbins = length(binedges) - 1
    bincenters = binedges[1:end-1] .+ 0.5*binwidth


    binheights = 1/binwidth*(h.weights .+ 1)/(length(data) + numbins + 1)
    sigma = (1/binwidth)*sqrt.( (1/(length(data) + numbins + 2))*(binwidth*binheights.*(1 .- binwidth*binheights)) )


    #fig = Figure(); display(fig)
    Axis(fig[2, 1], title="Distribution of density")
    

    barplot!(bincenters, binheights, width=binwidth,gap=0,label="MD simulation")
    errorbars!(bincenters, binheights, Vector{Float64}(vec(2*sigma)), color = :red)

    # define model for fit

    m(h,c) = c[1]*exp.(-h/c[2]) # fit parameter(s) have to be of type Array
    c0 = [10.0,10.0] # initial guess is important !!!



    fit_func = LsqFit.curve_fit(m,bincenters,binheights,c0)


    x_fit = LinRange(minimum(binedges), maximum(binedges),100)
    y_fit = m(x_fit,fit_func.param)

    plot!(x_fit,y_fit,label="LS-Fit")
    axislegend()

    @printf("kT/m from density histogram: %.3f \n", fit_func.param[2]*g)

    return nothing

end

function energy_timeseries(data) 

    Axis(fig[3, 1], title="Energy time series")
    plot!(data,label="Energy")
    axislegend()

    return nothing
end

function force(x,y,g)

    N = length(x)

    fx = zeros(N)
    fy = zeros(N)

    for i = 1:N-1
        for j = (i+1):N

            dx = x[i] - x[j]
            dy = y[i] - y[j]

            r_2 = dx^2 + dy^2

            LJ = 48/r_2*(1/r_2^6 - 1/(2*r_2^3))

            fx[i] += dx*LJ # how all the j affect i
            fy[i] += dy*LJ
            
    
            fx[j] += -dx*LJ # how i affects all the j
            fy[j] += -dy*LJ

        end
    end

    fy = fy .- g # add contribution of gravity to the force (m = 1 in reduced units), chosen in y-direction 

    return fx,fy
   
end

function energy(x,y,vx,vy,g)

    N = length(x)

    E_LJ = 0

    for i = 1:N-1
        for j = (i+1):N

            dx = x[i] - x[j]
            dy = y[i] - y[j]

            r_2 = dx^2 + dy^2

            E_LJ +=  4*(1/r_2^6 - 1/r_2^3)

        end
    end

    E_g = g*sum(y)
    E_kin = 0.5*sum(vx.^2 + vy.^2)

    return E_LJ + E_g + E_kin

end

function main() 

    # reduced units
    """
    k_b = 1.380649e-23 # J/K
    m = 6.63e-26 # kg
    epsilon = 120*k_b # J
    sigma = 3.4e-10 # m
    t = sigma*np.sqrt(m/epsilon) # characteristic time
    k_b_T = 300*k_b/epsilon # kT at room  temperature in reduced units
    g = 9.81*sigma*m/epsilon # in reduced units

    @printf("Time scale for argon in reduced units: t = %.3f", t)
    @printf("kT for argon in reduced units: kT = %.3f",k_b_T)
    @printf("Gravitational force in reuced units: g = %.3f",g)

    """

    # define system parameter
    N = Int64(40)
    L = Int64(20)
    x,y,vx,vy = initialize(N,L) # add temperature to initialization????

    g = 80

    # define simulation parameter
    n_meas = Int64(1e3)
    n_skip = Int64(1e1)
    n_loop = n_meas*n_skip
    n_therm = Int64(round(0.3*n_loop))
    animation = false

    tau = 0.0007

    # initialize measurements

    Y = zeros(n_meas,N) # for density calculation
    VX = zeros(n_meas,N)
    VY =  zeros(n_meas,N)
    E = zeros(n_meas)


    x,y,vx,vy = thermalization(x,y,vx,vy,n_therm,g,L,tau)
    Y,VX,VY,E = simulation_loop(x,y,vx,vy,n_loop,Y,VX,VY,E,g,L,n_skip,tau,animation)

    
    V = sqrt.(VX.^2 + VY.^2)
    V = V[:]

    fontsize_theme = Theme(fontsize = 25) # set plot parameter
    set_theme!(fontsize_theme)

    global fig = Figure(); display(fig)
    
    velocity_histogram(V)
    density_histogram(Y[:],g)
    energy_timeseries(E)


    return nothing

end


main()
