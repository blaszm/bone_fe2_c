function elmt01_r!(re::PseudoBlockArray{Float64,1}, elmt::CellIterator, cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, grid::Grid, mp::CorticalBone, msp::MacroSimulationParameters, ue::Array{Float64,1}, ve::Array{Float64,1}, ae::Array{Float64,1}, state::Array{MaterialState{Float64,Array{Float64,1}},1}, tr::Tensor{1,3,Float64,3}, t::Int, mq::MacroQuantities) where {dim}
    # ue = element displacement solution vector (at nodes)
    # ve = element velocity solution vector (at nodes)
    # ue speichert knotenweise erst alle u, dann alle phi, dann A
    
    n_basefuncs_u = getnbasefunctions(cellvalues_u) # Anzahl Ansatzfunktionen
    n_basefuncs_ϕ = getnbasefunctions(cellvalues_phi)
    n_basefuncs_A = getnbasefunctions(cellvalues_A)
    
    n_gp = getnquadpoints(cellvalues_u) # Anzahl GP
    
    u▄, ϕ▄ ,A▄ = 1,2,3 # Blockindex für gekoppelte Probleme
    reinit!(cellvalues_u, elmt)
    reinit!(cellvalues_phi, elmt)
    reinit!(cellvalues_A, elmt)
    
    # unpack simulation parameters
    ctan = zeros(3)
    ε_macro = zeros(6)
    E_macro = zeros(3)
    B_macro = zeros(3)
    
    ctan[:] = msp.ctan[:]
    γ = msp.γ
    ε_macro[:] = mq.ε_macro[:]
    E_macro[:] = mq.E_macro[:]
    B_macro[:] = mq.B_macro[:]
    
    
    # preallocate / reset helper variables
    Cᵉ = mp.Cᵉ 
    ϵₜ = mp.ϵₜ
    μⁱ = mp.μⁱ
    eₚ = mp.eₚ
    
    # mechanic
    ϵ = zeros(6)
    ϵp = zeros(6)
    σ = zeros(6)
    Bu = zeros(6,n_basefuncs_u)
    BuT = zeros(n_basefuncs_u,6)
    
    # electric
    E = zeros(3)
    D = zeros(3)
    Dp = zeros(3)
    Bgrad = zeros(3,n_basefuncs_ϕ)
    BgradT = zeros(n_basefuncs_ϕ,3)
    
    # magnetic
    B = zeros(3)
    H = zeros(3)
    N_A = zeros(3,n_basefuncs_A)
    N_AT = zeros(n_basefuncs_A,3)
    Bcurl = zeros(3,n_basefuncs_A)
    BcurlT = zeros(n_basefuncs_A,3)
    Bdiv = zeros(n_basefuncs_A)
    B_divdyad = zeros(n_basefuncs_A,n_basefuncs_A)
    div_Ψ = zeros(n_basefuncs_A)
    
    # residual vectors
    re_u = zeros(n_basefuncs_u)
    re_ϕ = zeros(n_basefuncs_ϕ)
    re_A = zeros(n_basefuncs_A)
    
    # main variables (nodes)
    u = zeros(n_basefuncs_u)
    ϕ = zeros(n_basefuncs_ϕ)
    A = zeros(n_basefuncs_A)
    up = zeros(n_basefuncs_u)
    ϕp = zeros(n_basefuncs_ϕ) 
    Ap = zeros(n_basefuncs_A)
    App = zeros(n_basefuncs_A)
    
    u[:] = ue[1:(n_basefuncs_u)]  
    ϕ[:] = ue[(n_basefuncs_u+1):(n_basefuncs_u+n_basefuncs_ϕ)] 
    A[:] = ue[(n_basefuncs_u+n_basefuncs_ϕ+1):(n_basefuncs_u+n_basefuncs_ϕ+n_basefuncs_A)]
    up[:] = ve[1:(n_basefuncs_u)]
    ϕp[:] = ve[(n_basefuncs_u+1):(n_basefuncs_u+n_basefuncs_ϕ)] 
    Ap[:] = ve[(n_basefuncs_u+n_basefuncs_ϕ+1):(n_basefuncs_u+n_basefuncs_ϕ+n_basefuncs_A)]
    App[:] = ae[(n_basefuncs_u+n_basefuncs_ϕ+1):(n_basefuncs_u+n_basefuncs_ϕ+n_basefuncs_A)]
    
    @inbounds for GPi in 1:n_gp # loop GPi
        
        # Bu operator matrix        
        Bu = voigtBu(cellvalues_u, GPi)
        BuT = Bu'
        
        # Bphi, Bcurl, NA, Bdiv 
        @inbounds for i in 1:n_basefuncs_A # basefuncs_u geht bis 24 wegen 3d
            N_A[:,i] = shape_value(cellvalues_A, GPi, i) 
            Bdiv[i] = shape_divergence(cellvalues_A, GPi, i)
            Bcurl[:,i] = shape_curl(cellvalues_A, GPi, i)
        end 
        
        @inbounds for i in 1:n_basefuncs_ϕ # bis 8
            Bgrad[:,i] = shape_gradient(cellvalues_phi, GPi, i) 
        end
        
        N_AT = N_A'
        BgradT = Bgrad'
        
        # volume part
        dΩ = getdetJdV(cellvalues_u, GPi)
        
        # strain + time derivative
        ε = voigtstrain(cellvalues_u, GPi, u)
        ε += ε_macro
        εp = voigtstrain(cellvalues_u, GPi, up)
        
        # electric field
        E = -1.0*(Bgrad * ϕ) - 1.0*(N_A * Ap)
        E += E_macro
        
        # magnetic flux density
        B = Bcurl * A
        B += B_macro
        
        # stress 
        σ = Cᵉ * ε - eₚ' * E
        
        # electric displacement field + time derivative
        D = ϵₜ * E + eₚ * ε
        Dp = ϵₜ * (-1.0*(Bgrad * ϕp) - 1.0*(N_A * App)) + eₚ * εp
        
        # magnetic field strength
        H = μⁱ * B
        
        # gauge function
        B_divdyad = Bdiv * Bdiv'
        div_Ψ = B_divdyad*A
        
        # save material state
        state[GPi].temp_σ = σ 
        state[GPi].temp_εⁱ = zeros(6)
        state[GPi].ε = ε
        state[GPi].E = E
        state[GPi].D = D
        state[GPi].Dp = Dp
        state[GPi].B = B
        state[GPi].H = H
        state[GPi].J = zeros(3)
        
        # Residual
        re_u -= (BuT * σ) * dΩ # altn deepcopy!
        re_ϕ -= (BgradT * D) * dΩ
        re_A -= (BcurlT * H - N_AT * Dp + γ * div_Ψ) * dΩ
        
    end # of loop GPi

    # Reassign to re with BlockIndex
    @inbounds for i in 1:(n_basefuncs_u)
        re[BlockIndex((u▄), (i))] = re_u[i]
    end
    
    @inbounds for i in 1:(n_basefuncs_A)
        re[BlockIndex((A▄), (i))] = re_A[i]
    end
    
    @inbounds for i in 1:n_basefuncs_ϕ
        re[BlockIndex((ϕ▄), (i))] = re_ϕ[i] 
    end
end

function elmt01_S!(Se::PseudoBlockArray{Float64,2}, elmt::CellIterator, cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, grid::Grid, mp::CorticalBone, msp::MacroSimulationParameters) where {dim}
    # ue = element displacement solution vector (at nodes)
    # ve = element velocity solution vector (at nodes)
    # ue speichert knotenweise erst alle u, dann alle phi, dann A
    
    n_basefuncs_u = getnbasefunctions(cellvalues_u) # Anzahl Ansatzfunktionen
    n_basefuncs_ϕ = getnbasefunctions(cellvalues_phi)
    n_basefuncs_A = getnbasefunctions(cellvalues_A)
    
    n_gp = getnquadpoints(cellvalues_u) # Anzahl GP
    
    u▄, ϕ▄ ,A▄ = 1,2,3 # Blockindex für gekoppelte Probleme
    reinit!(cellvalues_u, elmt)
    reinit!(cellvalues_phi, elmt)
    reinit!(cellvalues_A, elmt)
    
    # unpack simulation parameters
    ctan = zeros(3)
    ctan[:] = msp.ctan[:]
    γ = msp.γ
    
    # preallocate / reset helper variables
    Cᵉ = mp.Cᵉ 
    ϵₜ = mp.ϵₜ
    μⁱ = mp.μⁱ
    eₚ = mp.eₚ
    
    Bu = zeros(6,n_basefuncs_u)
    BuT = zeros(n_basefuncs_u,6)
    
    # electric
    Bgrad = zeros(3,n_basefuncs_ϕ)
    BgradT = zeros(n_basefuncs_ϕ,3)
    
    # magnetic
    N_A = zeros(3,n_basefuncs_A)
    N_AT = zeros(n_basefuncs_A,3)
    Bcurl = zeros(3,n_basefuncs_A)
    BcurlT = zeros(n_basefuncs_A,3)
    Bdiv = zeros(n_basefuncs_A)
    B_divdyad = zeros(n_basefuncs_A,n_basefuncs_A)
    
    # stiffness tensors
    Ke_uu = zeros(n_basefuncs_u,n_basefuncs_u)
    Ke_uϕ = zeros(n_basefuncs_u,n_basefuncs_ϕ)
    Ke_ϕu = zeros(n_basefuncs_ϕ,n_basefuncs_u)
    Ke_ϕϕ = zeros(n_basefuncs_ϕ,n_basefuncs_ϕ)
    Ke_AA = zeros(n_basefuncs_A,n_basefuncs_A)
    
    # damping tensors
    Ce_uA = zeros(n_basefuncs_u,n_basefuncs_A)
    Ce_ϕA = zeros(n_basefuncs_ϕ,n_basefuncs_A)
    Ce_Au = zeros(n_basefuncs_A,n_basefuncs_u)
    Ce_Aϕ = zeros(n_basefuncs_A,n_basefuncs_ϕ)
    
    # mass tensors
    Me_AA = zeros(n_basefuncs_A,n_basefuncs_A)
    
    @inbounds for GPi in 1:n_gp # loop GPi
        
        # Bu operator matrix
        Bu = voigtBu(cellvalues_u, GPi)
        BuT = Bu'
        
        # Bphi, Bcurl, NA, Bdiv 
        @inbounds for i in 1:n_basefuncs_A # basefuncs_u geht bis 24 wegen 3d
            N_A[:,i] = shape_value(cellvalues_A, GPi, i) 
            Bcurl[:,i] = shape_curl(cellvalues_A, GPi, i)
            Bdiv[i] = shape_divergence(cellvalues_A, GPi, i)
        end 
        
        @inbounds for i in 1:n_basefuncs_ϕ # bis 8
            Bgrad[:,i] = shape_gradient(cellvalues_phi, GPi, i) 
        end
        
        N_AT = N_A'
        BgradT = Bgrad'
        
        # gauge function
        B_divdyad = Bdiv * Bdiv'
        
        # volume part
        dΩ = getdetJdV(cellvalues_u, GPi)
        
        # stiffness tensors
        Ke_uu += (BuT * Cᵉ * Bu) * dΩ  
        Ke_uϕ += (BuT * eₚ' * Bgrad) * dΩ 
        Ke_ϕu += (BgradT * eₚ * Bu) * dΩ 
        Ke_ϕϕ += -1.0 * (BgradT * ϵₜ * Bgrad) * dΩ 
        Ke_AA += (BcurlT * μⁱ * Bcurl + γ * B_divdyad) * dΩ
        
        # damping tensors
        Ce_uA += (BuT * eₚ' * N_A) * dΩ
        Ce_ϕA += -1.0 * (BgradT * ϵₜ * N_A) * dΩ
        Ce_Au += -1.0 * (N_AT * eₚ * Bu) * dΩ
        Ce_Aϕ += (N_AT * ϵₜ * Bgrad) * dΩ
        
        # mass tensors
        Me_AA += (N_AT * ϵₜ * N_A) * dΩ
        
    end # of loop GPi

    # Reassign to Se with BlockIndex
    @inbounds for j in 1:(n_basefuncs_u)
        @inbounds for i in 1:(n_basefuncs_u)
            Se[BlockIndex((u▄,u▄), (i,j))] = ctan[1]*Ke_uu[i,j]
        end
    end
    
    @inbounds for j in 1:(n_basefuncs_A)
        @inbounds for i in 1:(n_basefuncs_u)
            Se[BlockIndex((u▄,A▄), (i,j))] = ctan[2]*Ce_uA[i,j]
        end
    end
    
    @inbounds for j in 1:(n_basefuncs_u)
        @inbounds for i in 1:(n_basefuncs_A)
            Se[BlockIndex((A▄,u▄), (i,j))] = ctan[2]*Ce_Au[i,j]
        end
    end
    
    @inbounds for j in 1:(n_basefuncs_A)
        @inbounds for i in 1:(n_basefuncs_A)
            Se[BlockIndex((A▄,A▄), (i,j))] = ctan[1]*Ke_AA[i,j] + ctan[3]*Me_AA[i,j]
        end
    end
    
    @inbounds for j in 1:(n_basefuncs_u)
        @inbounds for i in 1:n_basefuncs_ϕ
            Se[BlockIndex((ϕ▄,u▄), (i,j))] = ctan[1]*Ke_ϕu[i,j]
        end
    end
    
    @inbounds for j in 1:(n_basefuncs_A)
        @inbounds for i in 1:n_basefuncs_ϕ
            Se[BlockIndex((ϕ▄,A▄), (i,j))] = ctan[2]*Ce_ϕA[i,j]
        end
    end
    
    @inbounds for j in 1:n_basefuncs_ϕ
        @inbounds for i in 1:(n_basefuncs_u)    
            Se[BlockIndex((u▄,ϕ▄), (i,j))] = ctan[1]*Ke_uϕ[i,j] 
        end
        @inbounds for i in 1:(n_basefuncs_A) 
            Se[BlockIndex((A▄,ϕ▄), (i,j))] = ctan[2]*Ce_Aϕ[i,j] 
        end
        @inbounds for i in 1:n_basefuncs_ϕ
            Se[BlockIndex((ϕ▄,ϕ▄), (i,j))] = ctan[1]*Ke_ϕϕ[i,j]
        end
    end
end