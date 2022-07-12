# In PV:   (x,y,z,xy,yz,xz)

# returns strain vector in voigt notation for current element / GP
function voigtstrain(cellvalues_u::CellVectorValues{dim}, GPi::Int, ue::Array{Float64,1}) where {dim}

    ϵ_t = symmetric(function_gradient(cellvalues_u, GPi, ue)) # Tensor object with current strains
    
    ϵ = zeros(6)
    
    # Umbelegung für Voigt-Notation
    ϵ[1] = ϵ_t[1,1] # ϵ_x
    ϵ[2] = ϵ_t[2,2] # ϵ_y
    ϵ[3] = ϵ_t[3,3] # ϵ_z
    ϵ[4] = ϵ_t[1,2] + ϵ_t[2,1] # 2ϵ_xy 
    ϵ[5] = ϵ_t[2,3] + ϵ_t[3,2] # 2ϵ_yz
    ϵ[6] = ϵ_t[1,3] + ϵ_t[3,1] # 2ϵ_xz
    
    return ϵ
end

# returns Bu operator matrix in voigt notation for current element / GP
function voigtBu(cellvalues_u::CellVectorValues{dim}, GPi::Int) where {dim}
    n_basefuncs_u = getnbasefunctions(cellvalues_u) # geht bis 24 wegen 3d!
    n_basefuncs = Int(n_basefuncs_u/3.0) # 3.0 dim
    
    Bu = zeros(6,n_basefuncs_u)
    
    @inbounds for i in 1:n_basefuncs
        gNN = shape_gradient(cellvalues_u, GPi, 3*i-2)[1,:]
        
        Bu[1,1+(i-1)*3] = gNN[1] # eps_xx 
        Bu[2,2+(i-1)*3] = gNN[2] # eps_yy
        Bu[3,3+(i-1)*3] = gNN[3] # eps_zz
        
        Bu[4,1+(i-1)*3] = gNN[2] # 2*eps_xy  
        Bu[4,2+(i-1)*3] = gNN[1]
        Bu[5,2+(i-1)*3] = gNN[3] # 2*eps_yz
        Bu[5,3+(i-1)*3] = gNN[2]
        Bu[6,1+(i-1)*3] = gNN[3] # 2*eps_xz
        Bu[6,3+(i-1)*3] = gNN[1]
    end
    
    return Bu
end