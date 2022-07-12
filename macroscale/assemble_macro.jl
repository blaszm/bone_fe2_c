function doassemble_macro(S::SparseMatrixCSC, grid::Grid, dh::DofHandler, mp_b::CorticalBone, mp_m::BoneMarrow, sp::MicroProblem, mtm::MacroTangentModuli, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}, tr::Tensor{1,3,Float64,3}, t::Int, etype_macro::DataType, etype_micro::DataType)

    @timeit "parloop" result = pmap(i -> elmt_par!(i, grid, mp_b, mp_m, sp, mtm, u, v, a, tr, t, dh, etype_macro, etype_micro), 1:length(CellIterator(dh))) 
    
    r = zeros(ndofs(dh))
    assembler = start_assemble(S, r)
    for i in 1:length(CellIterator(dh))
        # unpack re, se, eldofs 
        se = result[i][1]
        re = result[i][2]
        eldofs = result[i][3]
        # assemble
        assemble!(assembler, eldofs, re, se)
    end
    
    # states
    for (el, state) in enumerate(states) 
        state_el = 0
        for i in 1:length(CellIterator(dh)) 
            if(el==result[i][5])
                state_el = i
                break
            end
        end
        for j in 1:length(state)
            state[j].σ[:] = result[state_el][4][j].σ[:]
            state[j].ε[:] = result[state_el][4][j].ε[:]
            state[j].E[:] = result[state_el][4][j].E[:]
            state[j].D[:] = result[state_el][4][j].D[:]
            state[j].Dp[:] = result[state_el][4][j].Dp[:]
            state[j].B[:] = result[state_el][4][j].B[:]
            state[j].H[:] = result[state_el][4][j].H[:]
            state[j].J[:] = result[state_el][4][j].J[:]
        end
    end 
    
    result = nothing
    
    return S, r
end;

function elmt_par!(i::Int, grid::Grid, mp_b::CorticalBone, mp_m::BoneMarrow, sp::MicroProblem, mtm::MacroTangentModuli, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, tr::Tensor{1,3,Float64,3}, t::Int, dh::DofHandler, etype_macro::DataType, etype_micro::DataType) where {dim}
    # i -> elmtno
    
    cellvalues_u, cellvalues_phi, cellvalues_A, facevalues_u = create_values(etype_macro)
    
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_phi)
    nA = getnbasefunctions(cellvalues_A)

    re = PseudoBlockArray(zeros(nu+nϕ+nA), [nu, nϕ, nA]) # local force vector
    se = PseudoBlockArray(zeros(nu+nϕ+nA,nu+nϕ+nA), [nu, nϕ, nA], [nu, nϕ, nA]) # local stiffness matrix
    fill!(se, 0)
    fill!(re, 0)
    
    eldofs = zeros(Int, ndofs_per_cell(dh))
    celldofs!(eldofs, dh, i)
    
    ue = zeros(ndofs_per_cell(dh))
    ve = zeros(ndofs_per_cell(dh))
    ae = zeros(ndofs_per_cell(dh))
    ue[:] = u[eldofs][:]
    ve[:] = v[eldofs][:]
    ae[:] = a[eldofs][:]
    
    n_gp = getnquadpoints(cellvalues_u)
    state = [MacroMaterialState() for _ in 1:n_gp]
    
    elmt_macro!(se, re, i, cellvalues_u, cellvalues_phi, cellvalues_A, facevalues_u, grid, mp_b, mp_m, sp, mtm, ue, ve, ae, state, tr, t, etype_micro)

    return se, re, eldofs, state, i
end