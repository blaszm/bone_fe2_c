# Simple grid creation
#=
function create_cook_grid(elemtype::DataType, nx::Int, ny::Int, nz::Int)
    # Ecken zuerst mit Koordinaten
    c1 = Vec{3}((-0.025, -0.015, -0.01)) # kleiner machen! RVE hat jetzt V = 0.001 m^3
    c2 = Vec{3}((0.025, 0.015, 0.01))
    grid = generate_grid(elemtype, (nx, ny, nz), c1, c2);
    addnodeset!(grid, "llcorner", x -> x[1] ≈ -0.025 && x[2] ≈ -0.015 && x[3] ≈ -0.01)
    return grid
end; 
=#

function create_bc_macro(dh::DofHandler, amplitude::Float64; modeltype::String="", option::String="default")
    dbc = ConstraintHandler(dh)

    if(modeltype == "cylinder")
        if(option == "default")
            add!(dbc, Dirichlet(:phi, getfaceset(dh.grid, "leftgrounding"), (x,t) -> zero(Vec{1}), [1]))
            add!(dbc, Dirichlet(:A, getfaceset(dh.grid, "leftgrounding"), (x,t) -> zero(Vec{3}), [1,2,3]))

            add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "leftsurf"), (x,t) -> zero(Vec{3}), [1,2,3]))
            add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "rightsurf"), (x,t) -> zero(Vec{3}), [1,2,3]))
            add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "middlesurf"), (x,t) -> amplitude*time_magnitude(t, n_timesteps), [1]))
        elseif(option == "torsion")
            add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "leftsurf"), (x,t) -> zero(Vec{3}), [1,2,3]))

            add!(dbc, Dirichlet(:phi, getfaceset(dh.grid, "leftgrounding"), (x,t) -> zero(Vec{1}), [1]))
            add!(dbc, Dirichlet(:A, getfaceset(dh.grid, "leftgrounding"), (x,t) -> zero(Vec{3}), [1,2,3]))

            add!(dbc, Dirichlet(:u, getnodeset(dh.grid, "torstop"), (x,t) -> amplitude*time_magnitude(t, n_timesteps), [2])) # 0.002 original
            add!(dbc, Dirichlet(:u, getnodeset(dh.grid, "torsbottom"), (x,t) -> -amplitude*time_magnitude(t, n_timesteps), [2]))  
        end
    elseif(modeltype == "femurbonescaled")
        #add!(dbc, Dirichlet(:phi, getnodeset(dh.grid, "leftsinglenode"), (x,t) -> zero(Vec{1}), [1]))
        #add!(dbc, Dirichlet(:A, getnodeset(dh.grid, "leftsinglenode"), (x,t) -> zero(Vec{3}), [1,2,3])) 
        
        add!(dbc, Dirichlet(:phi, getfaceset(dh.grid, "leftgrounding"), (x,t) -> zero(Vec{1}), [1]))
        add!(dbc, Dirichlet(:A, getfaceset(dh.grid, "leftgrounding"), (x,t) -> zero(Vec{3}), [1,2,3]))
        
        add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "leftsurf"), (x,t) -> zero(Vec{3}), [1,2,3]))
        add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "rightsurf"), (x,t) -> zero(Vec{3}), [1,2,3]))
        add!(dbc, Dirichlet(:u, getfaceset(dh.grid, "middlesurf"), (x,t) -> amplitude*time_magnitude(t, n_timesteps), [1]))
    end
    
    close!(dbc)
    t = 0.0 # t = time
    Ferrite.update!(dbc, t)
    return dbc
end # Kräfte seperat

function time_magnitude(time::Float64, maxtime::Int) # 0 -> 1 -> -1 -> 0 -cycle - maxtime-1 durch 4 teilbar!
    @assert mod(maxtime-1,4)==0
    t = Int(time)
    al = Int((maxtime-1)/4+1)
    bl = Int((maxtime-1)/2+1)
    cl = Int((maxtime-1)/4+1)
    a = range(0.0,1.0,length=al)
    b = range(1.0,-1.0,length=bl)
    c = range(-1.0,0.0,length=cl)
    if(t==0)
        f_m = 0.0
    elseif(t < al)
        f_m = a[t]
    elseif(t < (al+bl-1))
        f_m = b[t-al+1]
    else
        f_m = c[t-(al+bl-2)]
    end
    return f_m 
end