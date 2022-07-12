function doassemble_r(cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, grid::Grid, dh::DofHandler, mp_b::CorticalBone, mp_m::BoneMarrow, msp::MacroSimulationParameters, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}, tr::Tensor{1,3,Float64,3}, t::Int, mq::MacroQuantities) where {dim}

    r = zeros(ndofs(dh))
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_phi)
    nA = getnbasefunctions(cellvalues_A)
    re = PseudoBlockArray(zeros(nu+nϕ+nA), [nu, nϕ, nA]) # local force vector

    @inbounds for (elmt, state) in zip(CellIterator(dh), states)

        fill!(re, 0)
        eldofs = celldofs(elmt)
        ue = u[eldofs]
        ve = v[eldofs]
        ae = a[eldofs]
        
        if cellid(elmt) in getcellset(grid, "Material_1")
            elmt01_r!(re, elmt, cellvalues_u, cellvalues_phi, cellvalues_A, grid, mp_b, msp, ue, ve, ae, state, tr, t, mq)
        elseif cellid(elmt) in getcellset(grid, "Material_2")
            elmt02_r!(re, elmt, cellvalues_u, cellvalues_phi, cellvalues_A, grid, mp_m, msp, ue, ve, ae, state, tr, t, mq)
        else
            print("Error: Element without Material found")
        end
        
        assemble!(r, eldofs, re)
    end
    return r
end

function doassemble_S(cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim},
                    S::SparseMatrixCSC, grid::Grid, dh::DofHandler, mp_b::CorticalBone, mp_m::BoneMarrow, msp::MacroSimulationParameters) where {dim}
    
    assembler = start_assemble()
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_phi)
    nA = getnbasefunctions(cellvalues_A)
    Se = PseudoBlockArray(zeros(nu+nϕ+nA,nu+nϕ+nA), [nu, nϕ, nA], [nu, nϕ, nA]) # local stiffness matrix

    @inbounds for elmt in CellIterator(dh)
        fill!(Se, 0)
        eldofs = celldofs(elmt)
        if cellid(elmt) in getcellset(grid, "Material_1")
            elmt01_S!(Se, elmt, cellvalues_u, cellvalues_phi, cellvalues_A, grid, mp_b, msp)
        elseif cellid(elmt) in getcellset(grid, "Material_2")
            elmt02_S!(Se, elmt, cellvalues_u, cellvalues_phi, cellvalues_A, grid, mp_m, msp)
        else
            print("Error: Element without Material found")
        end
        
        assemble!(assembler, eldofs, Se)
    end
    S = end_assemble(assembler)
    return S
end