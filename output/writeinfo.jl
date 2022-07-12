# write file info.txt with simulation information
function writeinfo(mp_b::CorticalBone, mp_m::BoneMarrow, n_dim::Int, n_npe::Int, n_fpe::Int, n_dofs_n::Int, n_gp::Int, n_el::Int, Ω::Float64, ε_macro::Array{Float64,1}, E_macro::Array{Float64,1}, B_macro::Array{Float64,1}, σ_aver::Array{Float64,1}, D_aver::Array{Float64,1}, Dp_aver::Array{Float64,1}, H_aver::Array{Float64,1}, J_aver::Array{Float64,1}, el_mat1::Int, el_mat2::Int, globgpnumber::Int)
        work_dir = pwd()
        cd("Results")
    
        gpnumber = @sprintf("%5.5i", globgpnumber)
        info = string("info_", gpnumber, ".txt")
        io = open(info, "w") # ggbf anders noch falls close erst später kommt
        println(io, "Simulation information \n")
    
        println(io,"Material parameters \n")
    
        println(io,"Cortical bone \n")
        println(io,"E_b: $(@sprintf("%.4e", mp_b.E))")
        println(io,"nue_b: $(@sprintf("%.4e", mp_b.ν))")
        println(io,"eps_11_b: $(@sprintf("%.4e", mp_b.ϵ₁))")
        println(io,"eps_33_b: $(@sprintf("%.4e", mp_b.ϵ₃))")
        println(io,"mueinv_b: $(@sprintf("%.4e", mp_b.μᶜ))")
        println(io,"e_p14_b: $(@sprintf("%.4e", mp_b.eₚ14))")
        println(io, "\n")
    
        println(io, "Elastic stiffness tensor")    
        @inbounds for i in 1:length(mp_b.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.Cᵉ[i,1])) $(@sprintf("%.8e", mp_b.Cᵉ[i,2])) $(@sprintf("%.8e", mp_b.Cᵉ[i,3])) $(@sprintf("%.8e", mp_b.Cᵉ[i,4])) $(@sprintf("%.8e", mp_b.Cᵉ[i,5])) $(@sprintf("%.8e", mp_b.Cᵉ[i,6]))") 
        end
        println(io, "\n")
    
        println(io, "Dielectric tensor")  
        @inbounds for i in 1:length(mp_b.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.ϵₜ[i,1])) $(@sprintf("%.8e", mp_b.ϵₜ[i,2])) $(@sprintf("%.8e", mp_b.ϵₜ[i,3]))")
        end
        println(io, "\n")
    
        println(io, "Inverse magnetic tensor")
        @inbounds for i in 1:length(mp_b.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.μⁱ[i,1])) $(@sprintf("%.8e", mp_b.μⁱ[i,2])) $(@sprintf("%.8e", mp_b.μⁱ[i,3]))")
        end
        println(io, "\n")
    
        println(io, "Piezoelectric tensor")       
        @inbounds for i in 1:length(mp_b.eₚ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.eₚ[i,1])) $(@sprintf("%.8e", mp_b.eₚ[i,2])) $(@sprintf("%.8e", mp_b.eₚ[i,3])) $(@sprintf("%.8e", mp_b.eₚ[i,4])) $(@sprintf("%.8e", mp_b.eₚ[i,5])) $(@sprintf("%.8e", mp_b.eₚ[i,6]))") 
        end
    
        println(io, "\n")
        
        println(io,"Bone marrow \n")
        println(io,"E_m: $(@sprintf("%.4e", mp_m.E))")
        println(io,"nue_m: $(@sprintf("%.4e", mp_m.ν))")
        println(io,"eps_11_m: $(@sprintf("%.4e", mp_m.ϵ₁))")
        println(io,"eps_33_m: $(@sprintf("%.4e", mp_m.ϵ₃))")
        println(io,"mueinv_m: $(@sprintf("%.4e", mp_m.μᶜ))")
        println(io,"kappa_1_m: $(@sprintf("%.4e", mp_m.κ₁))")
        println(io,"mue_ve_m: $(@sprintf("%.4e", mp_m.μᵥ))")    
        println(io, "\n")
    
        println(io, "Elastic stiffness tensor")       
        @inbounds for i in 1:length(mp_m.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.Cᵉ[i,1])) $(@sprintf("%.8e", mp_m.Cᵉ[i,2])) $(@sprintf("%.8e", mp_m.Cᵉ[i,3])) $(@sprintf("%.8e", mp_m.Cᵉ[i,4])) $(@sprintf("%.8e", mp_m.Cᵉ[i,5])) $(@sprintf("%.8e", mp_m.Cᵉ[i,6]))") 
        end
        println(io, "\n")
    
        println(io, "Dielectric tensor")  
        @inbounds for i in 1:length(mp_m.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.ϵₜ[i,1])) $(@sprintf("%.8e", mp_m.ϵₜ[i,2])) $(@sprintf("%.8e", mp_m.ϵₜ[i,3]))")
        end    
        println(io, "\n")
    
        println(io, "Inverse magnetic tensor")
        @inbounds for i in 1:length(mp_m.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.μⁱ[i,1])) $(@sprintf("%.8e", mp_m.μⁱ[i,2])) $(@sprintf("%.8e", mp_m.μⁱ[i,3]))")
        end
        println(io, "\n")
    
        println(io, "Conductivity tensor")
        @inbounds for i in 1:length(mp_m.κ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.κ[i,1])) $(@sprintf("%.8e", mp_m.κ[i,2])) $(@sprintf("%.8e", mp_m.κ[i,3]))")
        end
    
        println(io, "\n")
    
        println(io, "Mesh parameters \n")
        #println(io, "Mesh resolution: $n")
        println(io, "Element type: ")
        println(io, "dim $n_dim, nodes per element $n_npe, faces per element $n_fpe")
        println(io, "DoF per node: $n_dofs_n")
        println(io, "Gaußpoints per element: $n_gp")
        println(io, "Total number of elements: $n_el")
        println(io, "Total Volume: $(@sprintf("%.4e", Ω))")
        println(io, "\n")
    
        println(io, "Number Elements of Material 1: $el_mat1")
        println(io, "Number Elements of Material 2: $el_mat2")
        println(io, "\n")
    
        # macro quantities
        println(io, "Macro strain:")
        println(io,"$(@sprintf("%.8e", ε_macro[1])) $(@sprintf("%.8e", ε_macro[2])) $(@sprintf("%.8e", ε_macro[3])) $(@sprintf("%.8e", ε_macro[4])) $(@sprintf("%.8e", ε_macro[5])) $(@sprintf("%.8e", ε_macro[6]))")
        println(io, "\n")
    
        println(io, "Macro electric field:")
        println(io,"$(@sprintf("%.8e", E_macro[1])) $(@sprintf("%.8e", E_macro[2])) $(@sprintf("%.8e", E_macro[3]))")
        println(io, "\n")
    
        println(io, "Macro magnetic flux density:")
        println(io,"$(@sprintf("%.8e", B_macro[1])) $(@sprintf("%.8e", B_macro[2])) $(@sprintf("%.8e", B_macro[3]))")
        println(io, "\n")
    
        # volume averages
        println(io, "Volume average stress:")
        println(io,"$(@sprintf("%.8e", σ_aver[1])) $(@sprintf("%.8e", σ_aver[2])) $(@sprintf("%.8e", σ_aver[3])) $(@sprintf("%.8e", σ_aver[4])) $(@sprintf("%.8e", σ_aver[5])) $(@sprintf("%.8e", σ_aver[6]))")
        println(io, "\n")
    
        println(io, "Volume average electric displacement:")
        println(io,"$(@sprintf("%.8e", D_aver[1])) $(@sprintf("%.8e", D_aver[2])) $(@sprintf("%.8e", D_aver[3]))")
        println(io, "\n")
    
        println(io, "Volume average electric displacement time derivative:")
        println(io,"$(@sprintf("%.8e", Dp_aver[1])) $(@sprintf("%.8e", Dp_aver[2])) $(@sprintf("%.8e", Dp_aver[3]))")
        println(io, "\n")
    
        println(io, "Volume average magnetic field strength:")
        println(io,"$(@sprintf("%.8e", H_aver[1])) $(@sprintf("%.8e", H_aver[2])) $(@sprintf("%.8e", H_aver[3]))")
        println(io, "\n")
    
        println(io, "Volume average electric current density:")
        println(io,"$(@sprintf("%.8e", J_aver[1])) $(@sprintf("%.8e", J_aver[2])) $(@sprintf("%.8e", J_aver[3]))")
        println(io, "\n")

        close(io)
        cd(work_dir)
    return
end

# write file info.txt with simulation information
function writeparameter(mp_b::CorticalBone, mp_m::BoneMarrow, n_dim::Int, n_npe::Int, n_fpe::Int, n_dofs_n::Int, n_gp::Int, n_el::Int, msp::MacroSimulationParameters, modeltype::String, option::String, etype_macro::DataType, etype_micro::DataType)
    
        work_dir = pwd()
        cd("Results")

        io = open("modelinfo.txt", "w") 
        println(io, "Simulation information \n")
    
        println(io,"Material parameters \n")
    
        println(io,"Cortical bone \n")
        println(io,"E_b: $(@sprintf("%.4e", mp_b.E))")
        println(io,"nue_b: $(@sprintf("%.4e", mp_b.ν))")
        println(io,"eps_11_b: $(@sprintf("%.4e", mp_b.ϵ₁))")
        println(io,"eps_33_b: $(@sprintf("%.4e", mp_b.ϵ₃))")
        println(io,"mueinv_b: $(@sprintf("%.4e", mp_b.μᶜ))")
        println(io,"e_p14_b: $(@sprintf("%.4e", mp_b.eₚ14))")
        println(io, "\n")
    
        println(io, "Elastic stiffness tensor")    
        @inbounds for i in 1:length(mp_b.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.Cᵉ[i,1])) $(@sprintf("%.8e", mp_b.Cᵉ[i,2])) $(@sprintf("%.8e", mp_b.Cᵉ[i,3])) $(@sprintf("%.8e", mp_b.Cᵉ[i,4])) $(@sprintf("%.8e", mp_b.Cᵉ[i,5])) $(@sprintf("%.8e", mp_b.Cᵉ[i,6]))") 
        end
        println(io, "\n")
    
        println(io, "Dielectric tensor")  
        @inbounds for i in 1:length(mp_b.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.ϵₜ[i,1])) $(@sprintf("%.8e", mp_b.ϵₜ[i,2])) $(@sprintf("%.8e", mp_b.ϵₜ[i,3]))")
        end
        println(io, "\n")
    
        println(io, "Inverse magnetic tensor")
        @inbounds for i in 1:length(mp_b.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.μⁱ[i,1])) $(@sprintf("%.8e", mp_b.μⁱ[i,2])) $(@sprintf("%.8e", mp_b.μⁱ[i,3]))")
        end
        println(io, "\n")
    
        println(io, "Piezoelectric tensor")       
        @inbounds for i in 1:length(mp_b.eₚ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.eₚ[i,1])) $(@sprintf("%.8e", mp_b.eₚ[i,2])) $(@sprintf("%.8e", mp_b.eₚ[i,3])) $(@sprintf("%.8e", mp_b.eₚ[i,4])) $(@sprintf("%.8e", mp_b.eₚ[i,5])) $(@sprintf("%.8e", mp_b.eₚ[i,6]))") 
        end
    
        println(io, "\n")
        
        println(io,"Bone marrow \n")
        println(io,"E_m: $(@sprintf("%.4e", mp_m.E))")
        println(io,"nue_m: $(@sprintf("%.4e", mp_m.ν))")
        println(io,"eps_11_m: $(@sprintf("%.4e", mp_m.ϵ₁))")
        println(io,"eps_33_m: $(@sprintf("%.4e", mp_m.ϵ₃))")
        println(io,"mueinv_m: $(@sprintf("%.4e", mp_m.μᶜ))")
        println(io,"kappa_1_m: $(@sprintf("%.4e", mp_m.κ₁))")
        println(io,"mue_ve_m: $(@sprintf("%.4e", mp_m.μᵥ))")    
        println(io, "\n")
    
        println(io, "Elastic stiffness tensor")       
        @inbounds for i in 1:length(mp_m.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.Cᵉ[i,1])) $(@sprintf("%.8e", mp_m.Cᵉ[i,2])) $(@sprintf("%.8e", mp_m.Cᵉ[i,3])) $(@sprintf("%.8e", mp_m.Cᵉ[i,4])) $(@sprintf("%.8e", mp_m.Cᵉ[i,5])) $(@sprintf("%.8e", mp_m.Cᵉ[i,6]))") 
        end
        println(io, "\n")
    
        println(io, "Dielectric tensor")  
        @inbounds for i in 1:length(mp_m.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.ϵₜ[i,1])) $(@sprintf("%.8e", mp_m.ϵₜ[i,2])) $(@sprintf("%.8e", mp_m.ϵₜ[i,3]))")
        end    
        println(io, "\n")
    
        println(io, "Inverse magnetic tensor")
        @inbounds for i in 1:length(mp_m.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.μⁱ[i,1])) $(@sprintf("%.8e", mp_m.μⁱ[i,2])) $(@sprintf("%.8e", mp_m.μⁱ[i,3]))")
        end
        println(io, "\n")
    
        println(io, "Conductivity tensor")
        @inbounds for i in 1:length(mp_m.κ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.κ[i,1])) $(@sprintf("%.8e", mp_m.κ[i,2])) $(@sprintf("%.8e", mp_m.κ[i,3]))")
        end
    
        println(io, "\n")
    
        # Macro mesh parameters
        println(io, "Macro mesh parameters \n")
        println(io, "Element type: ")
        println(io, "dim $n_dim, nodes per element $n_npe, faces per element $n_fpe")
        println(io, "DoF per node: $n_dofs_n")
        println(io, "Gaußpoints per element: $n_gp")
        println(io, "Total number of elements: $n_el")
        println(io, "\n")
    
        # Numerical parameters
        println(io, "Numerical parameters \n")
        
        println(io, "Rho_inf: $(@sprintf("%.4e", msp.ρ_∞))")
        println(io, "Delta t: $(@sprintf("%.4e", msp.Δₜ))")
        println(io, "Newton tolerance: $(@sprintf("%.4e", msp.NEWTON_TOL))")
        println(io, "Gauge penalization parameter gamma: $(@sprintf("%.4e", msp.γ))")
        println(io, "\n")
    
        # model
        println(io, "Model type: $modeltype")
        println(io, "option: $option")
        println(io, "Element type macro: $etype_macro")
        println(io, "Element type micro: $etype_micro")
    
        close(io)
        cd(work_dir)
    return
end