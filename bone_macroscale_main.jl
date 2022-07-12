const n_timesteps = 101
include("_include.jl")

function main()
    if(nprocs() > 1)
        for worker in workers()
            rmprocs(worker)
        end
    end

    if(nprocs() == 1)
        addprocs(40)
    end
    @everywhere include("_include.jl")
    
    solve_bone_macro("cylinder", "rve30pc", 2*1e-3; option="default", etype_macro = Hexahedron, etype_micro = Hexahedron); # call main with model, which should be used
    
    for worker in workers()
        rmprocs(worker)
    end
end

function initmaterialparameters(Δₜ::Float64)
    # cortical bone
    E = 22.e9 # Pa
    ν = 0.32 # -
    ϵ₁ = 8.85e-12 # F/m
    ϵ₃ = 8.85e-12 # F/m
    μᶜ = (1.257e-6)^(-1) # m/H, Inverse!
    eₚ14 = 3.e-3 # As/m^2, d*C 

    mp_b = CorticalBone(E, ν, ϵ₁, ϵ₃, μᶜ, eₚ14)

    # bone marrow
    E = 2.e9 # Pa
    ν = 0.3 # - 
    ϵ₁ = 8.85e-12 # F/m
    ϵ₃ = 8.85e-12 # F/m
    μᶜ = (1.257e-6)^(-1) # m/H, Inverse!
    κ₁ = 1.e4 # 5.e-3 # S/m <- im Vergl. zu H/Dp viel kleiner!
    μᵥ = 0.5*1e-9*Δₜ # s/Pa
    
    mp_m = BoneMarrow(E, ν, ϵ₁, ϵ₃, μᶜ, κ₁, μᵥ)
    return mp_b, mp_m
end

function solve_bone_macro(modeltype::String, micromodel::String, amplitude::Float64; option::String="default", etype_macro::DataType=Hexahedron, etype_micro::DataType=Hexahedron) # main
    reset_timer!()
    
    # micro mesh parameters
    a = 0.00032
    b = 0.00036
    n = 2
    
    # preperation export
    if(isdir("Results") == false) # erstellt Verzeichnis Results falls noch nicht vorhanden
        mkdir("Results")
    end

    msp = MacroSimulationParameters(0.5, 1e-4, 1e-8, 1.0) # ρ_∞, Δₜ, NEWTON_TOL, γ
    # unpack time integration parameters
    ρ_∞ = msp.ρ_∞
    α_f = msp.α_f
    αₘ = msp.αₘ
    γₐ = msp.γₐ
    Δₜ = msp.Δₜ
    NEWTON_TOL = msp.NEWTON_TOL

    # macromesh - grid, dofhandler, boundary condition
    n_dim = 3 # spacial dimension
    n_npe = 0 # number of nodes per element
    n_fpe = 0 # number of faces per element
    if(etype_macro == Hexahedron)
        n_npe = 8 
        n_fpe = 6 
    elseif(etype_macro == Tetrahedron)    
        n_npe = 4 
        n_fpe = 4 
    else
        println("Error: macro element type not supported")
    end
    elemtype = Cell{n_dim,n_npe,n_fpe} 
    n_macro = 3 # mesh resolution
    
    grid, cell_type, cellset, faceset = load_model(modeltype; option = option)
 
    dh = create_dofhandler(grid, etype_macro)
    dbc = create_bc_macro(dh, amplitude; modeltype = modeltype, option = option)

    # cellvalues
    cellvalues_u, cellvalues_phi, cellvalues_A, facevalues_u = create_values(etype_macro)
    
    # Pre-allocate solution vectors, etc.
    n_dofs = ndofs(dh)  # total number of dofs
    n_dofs_n = Int(ndofs_per_cell(dh)/n_npe) # DoF per node
    n_el = getncells(dh.grid) # total number of elements
    n_gp = getnquadpoints(cellvalues_u) # number of GP per element
    
    # current step t_n / start (t=0)
    dn = zeros(n_dofs) # displacement solution vector
    dpn = zeros(n_dofs) # velocity solution vector / time derivative
    vn = zeros(n_dofs) 
    vpn = zeros(n_dofs)

    # new step t_n+1
    dn1 = zeros(n_dofs)
    dpn1 = zeros(n_dofs)
    vn1 = zeros(n_dofs)
    vpn1 = zeros(n_dofs)
    
    # helper variables
    u_r = zeros(n_dofs) # residual displacements
    v_r = zeros(n_dofs) # residual velocities
    a_r = zeros(n_dofs) # residual accelarations
    
    # 3: vpn .= 0 # f(dn, dpn, F0) # ggbf Berechnung nötig, falls F =/= 0 z.B (Gl. 4)
    # 4: vn .= dpn # nur im ersten Schritt - dp/dpn können =/= 0 sein
    
    Δd = zeros(n_dofs) # solution increment
    r = zeros(n_dofs)  # residual
    S = create_sparsity_pattern(dh); # tangent stiffness matrix

    # apply!(dn, dbc)
    norm_r = 1.0
    norm_r0 = 1.0
    
    # create Material states
    states = [[MacroMaterialState() for _ in 1:n_gp] for _ in 1:getncells(grid)]
    
    # micro material (once)
    mp_b, mp_m = initmaterialparameters(Δₜ)
    
    # create micro problem
    sp = MicroProblem(micromodel, a, b, n, msp, mp_b, mp_m, etype_micro)
    
    # Vorbereitung HDF-files
    hdfname_i = "epsidata_i.h5"
    hdfname_ii = "epsidata_ii.h5"
    n_cells = length(sp.ctx.grid.cells) # n_cells = n_el_RVE
    preparehdf(hdfname_i, hdfname_ii, n_cells, n_el, n_gp)
        
    # macro tangent moduli
    mtm = calculatemacrotangent(sp, mp_b, mp_m, etype_micro)
    
    # Total Volume 
    Ω = 0.
    for (elmt, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues_u, cell) 
        for GPi in 1:n_gp 
            Ω += getdetJdV(cellvalues_u, GPi)
        end
    end 
    
    tr = Vec((0.0, 0.0, 0.0))
    
    # pvd file init 
    pvd = initpvd_macro()
    
    writeparameter(mp_b, mp_m, n_dim, n_npe, n_fpe, n_dofs_n, n_gp, n_el, msp, modeltype, option, etype_macro, etype_micro)
    
    #apply!(dn, dbc)
        
    # Newton-Raphson loop
    print("\n Starting Newton iterations:\n")
    
    for timestep in 1:n_timesteps
        t = timestep # aktueller Zeitschritt
        @inbounds dn1[:] = dn[:]
        
        Ferrite.update!(dbc, t)
        apply!(dn1, dbc)

        apply_zero!(dpn, dbc)
        apply_zero!(vn, dbc)
        apply_zero!(vpn, dbc)
        
        #tr_n = t > 1 ? Vec((0.01*time_magnitude(t-1), 0.0, 0.0)) : Vec((0.0, 0.0, 0.0))
        #tr_n1 = Vec((0.01*time_magnitude(t), 0.0, 0.0))
        
        newton_itr = -1
        print("\n Time step @time = $timestep:\n")

        while true; newton_itr += 1
            
            if newton_itr > 10
                error("Reached maximum Newton iterations, aborting")
                break
            end
        
            # !! vn1/vpn1 erfüllen RB nicht?!
            @inbounds vn1 = αₘ/(α_f*γₐ*Δₜ) * (dn1-dn) + (γₐ-αₘ)/(γₐ*α_f) * dpn + (α_f-1)/α_f * vn # (Gl. 18)
            @inbounds vpn1 =  αₘ/(α_f*γₐ^2*Δₜ^2)*(dn1-dn) - vn/(α_f*γₐ*Δₜ) + (γₐ-1)/γₐ * vpn + (γₐ-αₘ)/(α_f*γₐ^2*Δₜ) * dpn # (Gl. 19)

            @inbounds u_r = α_f*dn1+(1-α_f)*dn # (Gl.9) - dn+α_f
            @inbounds v_r = α_f*vn1+(1-α_f)*vn # (Gl.10) - vn+α_f
            @inbounds a_r = αₘ*vpn1+(1-αₘ)*vpn # (Gl.8) - vpn+αₘ
            #@inbounds Fnα_f = α_f*tr_n1 + (1-α_f)*tr_n # ggbf. (Gl. 11) falls Kräfte angreifen <--- !!
            
            # assembly and solve
            S, r = doassemble_macro(S, grid, dh, mp_b, mp_m, sp, mtm, u_r, v_r, a_r, states, tr, t, etype_macro, etype_micro);

            apply_zero!(S, r, dbc)            
            norm_r = norm(r)
            
            if(newton_itr==0 && norm_r > 0.0)
                norm_r0 = norm_r
            end
            
            print("Iteration: $newton_itr \trel. residual: $(@sprintf("%.4e", norm_r/norm_r0))\n")

            if (norm_r/norm_r0) < NEWTON_TOL
                break
            elseif norm_r < NEWTON_TOL
                print("Absolute norm chosen: $(@sprintf("%.4e", norm_r))\n")
                break
            end  
                        
            Δd,flag,err,iter,resvec = KrylovMethods.bicgstb(S, r; maxIter = 1000)
            
            apply_zero!(Δd, dbc)
            @inbounds dn1 .+= Δd

        end # of loop while NR-Iteration
                 
        # Compute new accelarations
        @inbounds dpn1 = (dn1 - dn)/(γₐ*Δₜ)+(γₐ-1)/γₐ * dpn # (Gl. 14)
    
        # 17 Update t_n+1 -> t_n
        @inbounds dn[:] = dn1[:]
        @inbounds dpn[:] = dpn1[:]
        @inbounds vn[:] = vn1[:]
        @inbounds vpn[:] = vpn1[:] 
        
        # export 
        writeoutput_macro(cellvalues_u, cellvalues_phi, dh, dbc, states, grid, t, dn, vn, vpn, pvd, r, etype_macro)

        updatehdf(hdfname_i, hdfname_ii)
        
    end # of loop timestep  
    
    writebc(dbc, grid, "bc_macro")
    
    savepvd(pvd)
    
    print_timer(title = "Test $n_timesteps Timesteps", linechars = :ascii)
end;

main();