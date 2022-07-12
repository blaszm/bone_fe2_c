function calculatemacrotangent(sp::MicroProblem, mp_b::CorticalBone, mp_m::BoneMarrow, etype_micro::DataType)
    ε_macro_p = zeros(6)
    E_macro_p = zeros(3)
    B_macro_p = zeros(3)
    
    σ_aver_p = zeros(6,6)
    D_aver_peps = zeros(3,6)
    D_aver_pE = zeros(3,3)
    H_aver_p = zeros(3,3)
    J_aver_p = zeros(3,3)
    
    C = zeros(6,6)
    ϵₜ = zeros(3,3)
    μⁱ = zeros(3,3)
    eₚ = zeros(3,6)
    κ = zeros(3,3)
    
    for i in 1:6 
        ε_macro_p .= 0.0
        E_macro_p .= 0.0
        B_macro_p .= 0.0
        ε_macro_p[i] += 1e-8
        mq = MacroQuantities(ε_macro_p, E_macro_p, B_macro_p, 1)
        σ_aver_p[:,i], D_aver_peps[:,i], Dp, H, J = solve_RVE(mq, sp, mp_b, mp_m, 1, etype_micro; overwrite=false)
        C[:,i] = (σ_aver_p[:,i])/1e-8
        eₚ[:,i] = (D_aver_peps[:,i])/1e-8
    end
    
    for i in 1:3
        ε_macro_p .= 0.0
        E_macro_p .= 0.0
        B_macro_p .= 0.0
        E_macro_p[i] += 1e-8
        mq = MacroQuantities(ε_macro_p, E_macro_p, B_macro_p, 1)
        σ, D_aver_pE[:,i], Dp, H, J_aver_p[:,i] = solve_RVE(mq, sp, mp_b, mp_m, 1, etype_micro; overwrite=false)
        ϵₜ[:,i] = (D_aver_pE[:,i])/1e-8
        κ[:,i] = (J_aver_p[:,i])/1e-8
    end
    
    for i in 1:3
        ε_macro_p .= 0.0
        E_macro_p .= 0.0
        B_macro_p .= 0.0
        B_macro_p[i] += 1e-8
        mq = MacroQuantities(ε_macro_p, E_macro_p, B_macro_p, 1)
        σ, D, Dp, H_aver_p[:,i], J = solve_RVE(mq, sp, mp_b, mp_m, 1, etype_micro; overwrite=false)
        μⁱ[:,i] = (H_aver_p[:,i])/1e-8        
    end
    
    # Make tensors symmetric 
    for i in 1:6
        for j in (i+1):6
            if(abs(C[j,i]) < 10.0)
                C[j,i] = 0
            end
            C[i,j] = C[j,i]
        end
    end
    
    for j in 1:6
        for i in 1:3
            if(abs(eₚ[i,j]) < 1e-12)
                eₚ[i,j] = 0.0
            end
        end
    end
    
    mtm = MacroTangentModuli(C, ϵₜ, μⁱ, eₚ, κ) # collect in mtm
    
    return mtm
end