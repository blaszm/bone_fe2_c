function solve_RVE(mq::MacroQuantities, sp::MicroProblem, mp_b::CorticalBone, mp_m::BoneMarrow, globgpnumber::Int, etype_micro::DataType; overwrite::Bool=true)
    #print("\n Starting RVE program:\n")
    
    # unpack parameters
    ε_macro = mq.ε_macro
    E_macro = mq.E_macro
    B_macro = mq.B_macro
    t_macro = mq.t

    a = sp.a
    b = sp.b
    n = sp.n
    ctx = sp.ctx
    bd = sp.bd
    S_red = sp.S_red
    
    ρ_∞ = sp.msp.ρ_∞
    α_f = sp.msp.α_f
    αₘ = sp.msp.αₘ
    γₐ = sp.msp.γₐ
    Δₜ = sp.msp.Δₜ  
    NEWTON_TOL = sp.msp.NEWTON_TOL
    nodedoflist = sp.nodedoflist
    
    # HDF file names
    hdfname_i = "epsidata_i.h5"
    hdfname_ii = "epsidata_ii.h5"
    dataname = "epsi_" * string(globgpnumber)

    n_timesteps = 1 
    
    # grid, dofhandler, boundary condition
    n_dim = 3 # spacial dimension
    n_npe = 0 # number of nodes per element
    n_fpe = 0 # number of faces per element
    if(etype_micro == Hexahedron)
        n_npe = 8 
        n_fpe = 6 
    elseif(etype_micro == Tetrahedron)    
        n_npe = 4 
        n_fpe = 4
    else
        println("Error: micro element type not supported")
    end
    elemtype = Cell{n_dim,n_npe,n_fpe}
    
    dh = create_dofhandler(ctx.grid, etype_micro)
    cellvalues_u, cellvalues_phi, cellvalues_A, facevalues_u = create_values(etype_micro)
    dbc = create_bc(dh)
    
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
    
    # 3
    # vpn .= 0 # f(dn, dpn, F0) # ggbf Berechnung nötig, falls F =/= 0 z.B (Gl. 4)
    # 4
    # vn .= dpn # nur im ersten Schritt - dp/dpn können =/= 0 sein
    
    r = zeros(n_dofs) # residual

    apply!(dn, dbc)
    
    r_red = ones(Int(nBCDofs(ctx, bd)*n_dofs_n))
    norm_r = 1.0
    norm_r0 = 1.0
    
    # creater helper variables for visualization of different sets
    material_number = zeros(getncells(ctx.grid))
    el_mat1 = 0 # Anzahl Elemente Material_1
    el_mat2 = 0 # Anzahl Elemente Material_2
    @inbounds for cell in CellIterator(dh)
        if cellid(cell) in getcellset(ctx.grid, "Material_1")
            material_number[cellid(cell)] = 1
            el_mat1 += 1
        elseif cellid(cell) in getcellset(ctx.grid, "Material_2")
            material_number[cellid(cell)] = 2
            el_mat2 += 1
        end
    end
    
    # create Material states
    n_gp = getnquadpoints(cellvalues_u) # für alle gleich
    states = [[MaterialState() for _ in 1:n_gp] for _ in 1:getncells(ctx.grid)] 
    
    # read εⁱ from HDF file and pass to states
    εⁱ_hdf = readhdf(hdfname_i, dataname)
    @inbounds for (elmt, cell) in enumerate(CellIterator(dh)) 
        state = states[elmt]
        @inbounds for GPi in 1:n_gp
            j = 8*(elmt-1) + GPi
            state[GPi].εⁱ[:] = εⁱ_hdf[:,j]
        end
    end
    
    # Total Volume 
    Ω = 0.0
    @inbounds for (elmt, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues_u, cell) 
        @inbounds for GPi in 1:n_gp 
            Ω += getdetJdV(cellvalues_u, GPi) 
        end
    end  
    
    # output variables
    σ_aver = zeros(6)
    D_aver = zeros(3)
    Dp_aver = zeros(3)
    H_aver = zeros(3)
    J_aver = zeros(3) 
    
    # traction vector
    tr = Vec{3}((0.0, 0.0, 0.0)) 
    
    # pvd file init 
    # pvd = initpvd(globgpnumber)

    # Solver
    #print("\n Starting Newton iterations:\n")
    
    #println(globgpnumber)
    
    for timestep in 1:n_timesteps # time loop
        
        t = timestep # aktueller Zeitschritt
        newton_itr = -1
        #print("\n Time step @time = $timestep:\n")
        @inbounds dn1[:] = dn[:]
        
        while true; newton_itr += 1
            
            if newton_itr > 10
                error("Reached maximum Newton iterations, aborting")
                break
            end
        
            @inbounds vn1 = αₘ/(α_f*γₐ*Δₜ) * (dn1-dn) + (γₐ-αₘ)/(γₐ*α_f) * dpn + (α_f-1)/α_f * vn # (Gl. 18)
            @inbounds vpn1 =  αₘ/(α_f*γₐ^2*Δₜ^2)*(dn1-dn) - vn/(α_f*γₐ*Δₜ) + (γₐ-1)/γₐ * vpn + (γₐ-αₘ)/(α_f*γₐ^2*Δₜ) * dpn # (Gl. 19)

            @inbounds u_r = α_f*dn1+(1-α_f)*dn # (Gl.9) - dn+α_f
            @inbounds v_r = α_f*vn1+(1-α_f)*vn # (Gl.10) - vn+α_f
            @inbounds a_r = αₘ*vpn1+(1-αₘ)*vpn # (Gl.8) - vpn+αₘ
            # ggbf. (Gl. 11) falls Kräfte angreifen
            
            # assembly and solve
            r = doassemble_r(cellvalues_u, cellvalues_phi, cellvalues_A, ctx.grid, dh, mp_b, mp_m, sp.msp, u_r, v_r, a_r, states, tr, t, mq); 

            apply_zero!(r, dbc)
            
            Δd_red = zeros(Int(nBCDofs(ctx, bd)*n_dofs_n))

            r_red = applyPBC_r(ctx, bd, r, nodedoflist, n_dofs_n)
            
            norm_r = norm(r_red)
            
            if(newton_itr==0 && norm_r > 0.0)
                norm_r0 = norm_r
            end
            
            #print("Iteration: $newton_itr \trel. residual: $(@sprintf("%.4e", norm_r/norm_r0))\n")

            if (norm_r/norm_r0) < NEWTON_TOL
                break
            elseif norm_r < NEWTON_TOL
                #print("Absolute norm chosen: $(@sprintf("%.4e", norm_r))\n")
                break
            end  

            #@timeit "krylovmethods bicgstb" 
            Δd_red,flag,err,iter,resvec = KrylovMethods.bicgstb(S_red, r_red; maxIter = 250)

            Δd = undoPBC(ctx, bd, n_dofs, n_dofs_n, Δd_red, nodedoflist)
            apply_zero!(Δd, dbc)
            @inbounds dn1 .+= Δd

        end # of loop while NR-Iteration
            
        # Update all the material states after we have reached equilibrium
        for cell_states in states
            foreach(update_state!, cell_states) # including update inelastic strain (internal variable)
        end
        
        # Compute new accelarations
        @inbounds dpn1 = (dn1 - dn)/(γₐ*Δₜ)+(γₐ-1)/γₐ * dpn # (Gl. 14)
    
        # 17 Update t_n+1 -> t_n
        @inbounds dn[:] = dn1[:]
        @inbounds dpn[:] = dpn1[:]
        @inbounds vn[:] = vn1[:]
        @inbounds vpn[:] = vpn1[:] 
        
        # volume averages
        σ_aver, D_aver, Dp_aver, H_aver, J_aver = volaver(cellvalues_A, dh, n_gp, Ω, states)
    
        # export
        #if(globgpnumber == 1)
        #    @timeit "write output" writeoutput(cellvalues_u, cellvalues_phi, dh, dbc, states, ctx.grid, t_macro, dn, vn, vpn, material_number, globgpnumber, pvd, etype_micro)
        #end
        
    end # of loop timestep  
        
    # write boundary_conditions.vtu
    # writebc(dbc, ctx.grid, "bc_micro")

    # write info.txt
    #if(globgpnumber == 1)
    #    writeinfo(mp_b, mp_m, n_dim, n_npe, n_fpe, n_dofs_n, n_gp, n_el, Ω, ε_macro, E_macro, B_macro, σ_aver, D_aver, Dp_aver, H_aver, J_aver, el_mat1, el_mat2, globgpnumber)
    #end
    
    # collect and save new εⁱ in HDF file
    if(overwrite)
        @inbounds for (elmt, cell) in enumerate(CellIterator(dh)) 
            state = states[elmt]
            @inbounds for GPi in 1:n_gp
                j = 8*(elmt-1) + GPi
                εⁱ_hdf[:,j] = state[GPi].εⁱ[:]
            end
        end
        overwritehdf(hdfname_ii, dataname, εⁱ_hdf)
    end
    
    #savepvd(pvd)
    
    return σ_aver, D_aver, Dp_aver, H_aver, J_aver
end