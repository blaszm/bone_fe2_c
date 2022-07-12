function collect_stresses(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}, i::Int) where {dim,T} 

    # i = stress component, order: x,y,z,xy,yz,xz - same as PV default
    sigma_out = [Vec{1,T}[] for _ in 1:getncells(dh.grid)] # Array(Array(Tensor))
    sigma_i = [0.0] # Array
    
    @inbounds for (el, cell_states) in enumerate(states)
        sigma_cell = sigma_out[el]
        
        @inbounds for state in cell_states
            sigma_i[1] = state.σ[i]
            sigma_t = reinterpret(Vec{1, T}, vec(sigma_i)) # Umwandlung Array -> Tensor
            push!(sigma_cell, sigma_t[1]) # Anhängen des Tensors an Array(Array())
        end
    end
    
    return sigma_out
end

function collect_D(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T} 
    D_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    D_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        D_cell = D_out[el]
        
        @inbounds for state in cell_states
            D_i[:] = state.D[:]
            D_t = reinterpret(Vec{3, T}, vec(D_i[:]))
            push!(D_cell, D_t[1])
        end
    end
    
    return D_out
end

function collect_Dp(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    Dp_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    Dp_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        Dp_cell = Dp_out[el]
        
        @inbounds for state in cell_states
            Dp_i[:] = state.Dp[:]
            Dp_t = reinterpret(Vec{3, T}, vec(Dp_i[:]))
            push!(Dp_cell, Dp_t[1])
        end
    end
    
    return Dp_out
end

function collect_E(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    E_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    E_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        E_cell = E_out[el]
        
        @inbounds for state in cell_states
            E_i[:] = state.E[:]
            E_t = reinterpret(Vec{3, T}, vec(E_i[:]))
            push!(E_cell, E_t[1])
        end
    end
    
    return E_out
end

function collect_J(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    J_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    J_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        J_cell = J_out[el]
        
        @inbounds for state in cell_states
            J_i[:] = state.J[:]
            J_t = reinterpret(Vec{3, T}, vec(J_i[:]))
            push!(J_cell, J_t[1])
        end
    end
    
    return J_out
end

function collect_B(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    B_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    B_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        B_cell = B_out[el]
        
        @inbounds for state in cell_states
            B_i[:] = state.B[:]
            B_t = reinterpret(Vec{3, T}, vec(B_i[:]))
            push!(B_cell, B_t[1])
        end
    end
    
    return B_out
end

function collect_H(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    H_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    H_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        H_cell = H_out[el]
        
        @inbounds for state in cell_states
            H_i[:] = state.H[:]
            H_t = reinterpret(Vec{3, T}, vec(H_i[:]))
            push!(H_cell, H_t[1])
        end
    end
    
    return H_out
end

function collect_strain(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}, i::Int) where {dim,T}

    # i = strain component, order: x,y,z,xy,yz,xz - same as PV default
    eps_out = [Vec{1,T}[] for _ in 1:getncells(dh.grid)] # Array(Array(Tensor))
    eps_i = [0.0] # Array
    
    @inbounds for (el, cell_states) in enumerate(states)
        eps_cell = eps_out[el]
        
        @inbounds for state in cell_states
            eps_i[1] = state.ε[i]
            eps_t = reinterpret(Vec{1, T}, vec(eps_i)) # Umwandlung Array -> Tensor
            push!(eps_cell, eps_t[1]) # Anhängen des Tensors an Array(Array())
        end
    end
    
    return eps_out
end

function collect_inelstrain(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}, i::Int) where {dim,T} 

    # i = strain component, order: x,y,z,xy,yz,xz - same as PV default
    epsi_out = [Vec{1,T}[] for _ in 1:getncells(dh.grid)] # Array(Array(Tensor))
    epsi_i = [0.0] # Array
    
    @inbounds for (el, cell_states) in enumerate(states)
        epsi_cell = epsi_out[el]
        
        @inbounds for state in cell_states
            epsi_i[1] = state.εⁱ[i]
            epsi_t = reinterpret(Vec{1, T}, vec(epsi_i)) # Umwandlung Array -> Tensor
            push!(epsi_cell, epsi_t[1]) # Anhängen des Tensors an Array(Array())
        end
    end
    
    return epsi_out
end

# --- overloaded methods for MacroMaterialState
function collect_stresses(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}, i::Int) where {dim,T} 

    # i = stress component, order: x,y,z,xy,yz,xz - same as PV default
    sigma_out = [Vec{1,T}[] for _ in 1:getncells(dh.grid)] # Array(Array(Tensor))
    sigma_i = [0.0] # Array
    
    @inbounds for (el, cell_states) in enumerate(states)
        sigma_cell = sigma_out[el]
        
        @inbounds for state in cell_states
            sigma_i[1] = state.σ[i]
            sigma_t = reinterpret(Vec{1, T}, vec(sigma_i)) # Umwandlung Array -> Tensor
            push!(sigma_cell, sigma_t[1]) # Anhängen des Tensors an Array(Array())
        end
    end
    
    return sigma_out
end

function collect_D(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T} 
    D_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    D_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        D_cell = D_out[el]
        
        @inbounds for state in cell_states
            D_i[:] = state.D[:]
            D_t = reinterpret(Vec{3, T}, vec(D_i[:]))
            push!(D_cell, D_t[1])
        end
    end
    
    return D_out
end

function collect_Dp(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    Dp_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    Dp_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        Dp_cell = Dp_out[el]
        
        @inbounds for state in cell_states
            Dp_i[:] = state.Dp[:]
            Dp_t = reinterpret(Vec{3, T}, vec(Dp_i[:]))
            push!(Dp_cell, Dp_t[1])
        end
    end
    
    return Dp_out
end

function collect_E(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    E_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    E_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        E_cell = E_out[el]
        
        @inbounds for state in cell_states
            E_i[:] = state.E[:]
            E_t = reinterpret(Vec{3, T}, vec(E_i[:]))
            push!(E_cell, E_t[1])
        end
    end
    
    return E_out
end

function collect_J(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    J_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    J_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        J_cell = J_out[el]
        
        @inbounds for state in cell_states
            J_i[:] = state.J[:]
            J_t = reinterpret(Vec{3, T}, vec(J_i[:]))
            push!(J_cell, J_t[1])
        end
    end
    
    return J_out
end

function collect_B(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    B_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    B_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        B_cell = B_out[el]
        
        @inbounds for state in cell_states
            B_i[:] = state.B[:]
            B_t = reinterpret(Vec{3, T}, vec(B_i[:]))
            push!(B_cell, B_t[1])
        end
    end
    
    return B_out
end

function collect_H(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}) where {dim,T}
    H_out = [Vec{3,T}[] for _ in 1:getncells(dh.grid)]
    H_i = [0.0, 0.0, 0.0]
    
    @inbounds for (el, cell_states) in enumerate(states)
        H_cell = H_out[el]
        
        @inbounds for state in cell_states
            H_i[:] = state.H[:]
            H_t = reinterpret(Vec{3, T}, vec(H_i[:]))
            push!(H_cell, H_t[1])
        end
    end
    
    return H_out
end

function collect_strain(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}, i::Int) where {dim,T}

    # i = strain component, order: x,y,z,xy,yz,xz - same as PV default
    eps_out = [Vec{1,T}[] for _ in 1:getncells(dh.grid)] # Array(Array(Tensor))
    eps_i = [0.0] # Array
    
    @inbounds for (el, cell_states) in enumerate(states)
        eps_cell = eps_out[el]
        
        @inbounds for state in cell_states
            eps_i[1] = state.ε[i]
            eps_t = reinterpret(Vec{1, T}, vec(eps_i)) # Umwandlung Array -> Tensor
            push!(eps_cell, eps_t[1]) # Anhängen des Tensors an Array(Array())
        end
    end
    
    return eps_out
end

function collect_inelstrain(cellvalues::CellVectorValues{dim,T}, dh::DofHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}, i::Int) where {dim,T} 

    # i = strain component, order: x,y,z,xy,yz,xz - same as PV default
    epsi_out = [Vec{1,T}[] for _ in 1:getncells(dh.grid)] # Array(Array(Tensor))
    epsi_i = [0.0] # Array
    
    @inbounds for (el, cell_states) in enumerate(states)
        epsi_cell = epsi_out[el]
        
        @inbounds for state in cell_states
            epsi_i[1] = state.εⁱ[i]
            epsi_t = reinterpret(Vec{1, T}, vec(epsi_i)) # Umwandlung Array -> Tensor
            push!(epsi_cell, epsi_t[1]) # Anhängen des Tensors an Array(Array())
        end
    end
    
    return epsi_out
end