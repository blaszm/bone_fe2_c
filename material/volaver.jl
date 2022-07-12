# calculates all necessary volume averages
function volaver(cellvalues_A::CellVectorValues{dim}, dh::DofHandler, n_gp::Int, Ω::Float64, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}) where {dim}

    σ_aver = zeros(6)
    D_aver = zeros(3)
    Dp_aver = zeros(3)
    H_aver = zeros(3)
    J_aver = zeros(3)
    
    @inbounds for (elmt, cell) in enumerate(CellIterator(dh)) 
        reinit!(cellvalues_A, cell)
        state = states[elmt]
        @inbounds for GPi in 1:n_gp 
            dΩ =  getdetJdV(cellvalues_A, GPi)
            
            σ_aver += state[GPi].σ * dΩ
            D_aver += state[GPi].D * dΩ
            Dp_aver += state[GPi].Dp * dΩ
            H_aver += state[GPi].H * dΩ
            J_aver += state[GPi].J * dΩ
            
        end
    end

    σ_aver /= Ω
    D_aver /= Ω
    Dp_aver /= Ω
    H_aver /= Ω
    J_aver /= Ω
    
    return σ_aver, D_aver, Dp_aver, H_aver, J_aver
end