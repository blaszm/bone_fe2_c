function elmt_macro!(Se::PseudoBlockArray{Float64,2}, re::PseudoBlockArray{Float64,1}, elmtno::Int, cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, facevalues_u::FaceVectorValues{dim}, grid::Grid, mp_b::CorticalBone, mp_m::BoneMarrow, sp::MicroProblem, mtm::MacroTangentModuli, ue::Array{Float64,1}, ve::Array{Float64,1}, ae::Array{Float64,1}, state::Array{MacroMaterialState{Float64,Array{Float64,1}},1}, tr::Tensor{1,3,Float64,3}, t::Int, etype_micro::DataType) where {dim} 
    
    # ue = element displacement solution vector (at nodes)
    # ve = element velocity solution vector (at nodes)
    # ae = element accelaration solution vector (at nodes)
    # ue speichert knotenweise erst alle u, dann alle phi, dann A
    
    n_basefuncs_u = getnbasefunctions(cellvalues_u) # Anzahl Ansatzfunktionen
    n_basefuncs_ϕ = getnbasefunctions(cellvalues_phi)
    n_basefuncs_A = getnbasefunctions(cellvalues_A)
    
    n_gp = getnquadpoints(cellvalues_u) # Anzahl GP
    
    coordinates = [zero(Vec{3}) for i in 1:length(grid.cells[1].nodes)]
    
    nodeids = grid.cells[elmtno].nodes
    for j in 1:length(coordinates)
        coordinates[j] = grid.nodes[nodeids[j]].x
    end 
    
    u▄, ϕ▄ ,A▄ = 1,2,3 # Blockindex für gekoppelte Probleme
    reinit!(cellvalues_u, coordinates) # <- par Problem?
    reinit!(cellvalues_phi, coordinates)
    reinit!(cellvalues_A, coordinates)
    
    # unpack simulation parameters
    ctan = zeros(3)
    ctan[:] = sp.msp.ctan[:]
    γ = sp.msp.γ
    
    # macro tangent moduli
    C = mtm.C
    ϵₜ = mtm.ϵₜ
    μⁱ = mtm.μⁱ
    eₚ = mtm.eₚ
    κ = mtm.κ
    
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
    J = zeros(3)
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
    
    # stiffness tensors
    Ke_uu = zeros(n_basefuncs_u,n_basefuncs_u)
    Ke_uϕ = zeros(n_basefuncs_u,n_basefuncs_ϕ)
    Ke_ϕu = zeros(n_basefuncs_ϕ,n_basefuncs_u)
    Ke_ϕϕ = zeros(n_basefuncs_ϕ,n_basefuncs_ϕ)
    Ke_Aϕ = zeros(n_basefuncs_A,n_basefuncs_ϕ)
    Ke_AA = zeros(n_basefuncs_A,n_basefuncs_A)
    
    # damping tensors
    Ce_uA = zeros(n_basefuncs_u,n_basefuncs_A)
    Ce_ϕA = zeros(n_basefuncs_ϕ,n_basefuncs_A)
    Ce_Au = zeros(n_basefuncs_A,n_basefuncs_u)
    Ce_Aϕ = zeros(n_basefuncs_A,n_basefuncs_ϕ)
    Ce_AA = zeros(n_basefuncs_A,n_basefuncs_A)
    
    # mass tensors
    Me_AA = zeros(n_basefuncs_A,n_basefuncs_A)
    
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
    
    @inbounds for GPi in 1:getnquadpoints(cellvalues_u) # loop GPi
        
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
        
        # gauge function
        B_divdyad = Bdiv * Bdiv'
        div_Ψ = B_divdyad*A
        
        # volume part
        dΩ = getdetJdV(cellvalues_u, GPi)        
        
        # stiffness tensors
        Ke_uu += (BuT * C * Bu) * dΩ  
        Ke_uϕ += (BuT * eₚ' * Bgrad) * dΩ 
        Ke_ϕu += (BgradT * eₚ * Bu) * dΩ 
        Ke_ϕϕ += -1.0 * (BgradT * ϵₜ * Bgrad) * dΩ 
        Ke_AA += (BcurlT * μⁱ * Bcurl + γ * B_divdyad) * dΩ
        Ke_Aϕ += (N_AT * κ * Bgrad) * dΩ
        
        # damping tensors
        Ce_uA += (BuT * eₚ' * N_A) * dΩ
        Ce_ϕA += -1.0 * (BgradT * ϵₜ * N_A) * dΩ
        Ce_Au += -1.0 * (N_AT * eₚ * Bu) * dΩ
        Ce_Aϕ += (N_AT * ϵₜ * Bgrad) * dΩ
        Ce_AA += (N_AT * κ * N_A) * dΩ
        
        # mass tensors
        Me_AA += (N_AT * ϵₜ * N_A) * dΩ 
        
        # strain + time derivative
        ε = voigtstrain(cellvalues_u, GPi, u)
        εp = voigtstrain(cellvalues_u, GPi, up)
        
        # electric field
        E = -1.0*(Bgrad * ϕ) - 1.0*(N_A * Ap)
        
        # magnetic flux density
        B = Bcurl * A      
        
        # --- micro scale ---
        # material -> solve RVE
        mq = MacroQuantities(ε, E, B, t)
        globgpnumber = Int(getnquadpoints(cellvalues_u)*(elmtno-1) + GPi)
        #print("Starting calculation of RVE no: ", globgpnumber, "\n")
        @timeit "solveRVE" σ, D, Dp, H, J = solve_RVE(mq, sp, mp_b, mp_m, globgpnumber, etype_micro) 
        #print("Received results from RVE no: ", globgpnumber, "\n")
        
        # save material state
        state[GPi].σ = σ 
        state[GPi].ε = ε
        state[GPi].E = E
        state[GPi].D = D
        state[GPi].Dp = Dp
        state[GPi].B = B
        state[GPi].H = H
        state[GPi].J = J
        
        # Residual
        re_u -= (BuT * σ) * dΩ # altn deepcopy!
        re_ϕ -= (BgradT * D) * dΩ
        re_A -= (BcurlT * H - N_AT * Dp - N_AT * J + γ * div_Ψ) * dΩ
        
    end # of loop GPi
    
    #=
    # Add traction as a negative contribution to the element residual `re`:
    for face in 1:nfaces(elmt)
        if onboundary(elmt, face) && (cellid(elmt), face) ∈ getfaceset(grid, "right") # hier verschiedene Möglichkeiten angeben
            reinit!(facevalues_u, elmt, face)
            for GPi in 1:getnquadpoints(facevalues_u)
                dΓ = getdetJdV(facevalues_u, GPi)
                for i in 1:n_basefuncs_u
                    δu = shape_value(facevalues_u, GPi, i)
                    re_u[i] += (δu ⋅ tr) * dΓ
                end
            end
        end
    end
    =#
    
    # Reassign to re / Se with BlockIndex
    @inbounds for j in 1:(n_basefuncs_u)
        re[BlockIndex((u▄), (j))] = re_u[j]
        @inbounds for i in 1:(n_basefuncs_u)
            Se[BlockIndex((u▄,u▄), (i,j))] = ctan[1]*Ke_uu[i,j]
        end
        @inbounds for i in 1:(n_basefuncs_A)
            Se[BlockIndex((A▄,u▄), (i,j))] = ctan[2]*Ce_Au[i,j]
        end
    end
    
    @inbounds for j in 1:(n_basefuncs_A)
        re[BlockIndex((A▄), (j))] = re_A[j]
        @inbounds for i in 1:(n_basefuncs_A)
            Se[BlockIndex((A▄,A▄), (i,j))] = ctan[1]*Ke_AA[i,j] + ctan[2]*Ce_AA[i,j] + ctan[3]*Me_AA[i,j]
        end
        @inbounds for i in 1:(n_basefuncs_u)
            Se[BlockIndex((u▄,A▄), (i,j))] = ctan[2]*Ce_uA[i,j]
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
        re[BlockIndex((ϕ▄), (j))] = re_ϕ[j]
        @inbounds for i in 1:(n_basefuncs_u)
            Se[BlockIndex((u▄,ϕ▄), (i,j))] = ctan[1]*Ke_uϕ[i,j] 
        end
        @inbounds for i in 1:(n_basefuncs_A)
            Se[BlockIndex((A▄,ϕ▄), (i,j))] = ctan[1]*Ke_Aϕ[i,j] + ctan[2]*Ce_Aϕ[i,j] 
        end
        @inbounds for i in 1:n_basefuncs_ϕ
            Se[BlockIndex((ϕ▄,ϕ▄), (i,j))] = ctan[1]*Ke_ϕϕ[i,j]
        end
    end   
    
end