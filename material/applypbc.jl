# (c) Mischa Blaszczyk 2020 - Ruhr-Universität Bochum
#  Email: mischa.blaszczyk@rub.de

function nodedofs(dh::DofHandler, grid::Grid) # nodedofs(i, :) = [nodenumber, fielddof], where i = dof in dh
    fielddims = zeros(Int64, Ferrite.nfields(dh)) 
    offsets = zeros(Int64, Ferrite.nfields(dh)) 
    
    @inbounds for (i, fielddim) in enumerate(dh.field_dims)
        fielddims[i] = fielddim
    end
    
    n_dofs = ndofs(dh)  # total number of dofs
    nnodes = length(grid.nodes) 
    nodedoflist = zeros(Int, n_dofs, 2)
    
    @inbounds for (i,fieldname) in enumerate(dh.field_names)
        offsets[i] = Ferrite.field_offset(dh, fieldname)
    end
    
    @inbounds for (i,node) in enumerate(grid.nodes)
        @inbounds for cell in CellIterator(dh)
            if(i in cell.nodes)
                cell_node = 0
                cell_dofs = celldofs(cell)
                @inbounds for j in 1:length(cell.nodes)
                    if(i == cell.nodes[j])
                        cell_node = j
                        break
                    end
                end
                k = 1
                @inbounds for j in 1:Ferrite.nfields(dh)
                    offset = offsets[j]
                    fielddim = fielddims[j]
                    if(j>1)
                        k += fielddims[j-1]
                    end
                    nodedoflist[cell_dofs[Int((cell_node-1)*fielddim + 1 + offset):Int((cell_node-1)*fielddim + offset + fielddim)], 1] .= i
                    nodedoflist[cell_dofs[Int((cell_node-1)*fielddim + 1 + offset):Int((cell_node-1)*fielddim + offset + fielddim)], 2] = Int(k):Int(k+fielddim-1) 
                end
                break
            end
        end
    end

    return nodedoflist
end

function applyPBC_S(ctx::CoherentStructures.GridContext, bdata::CoherentStructures.BoundaryData, S::SparseMatrixCSC, nodedoflist::Array{Int64,2}, ndof::Int) # gridcontext/boundary data/stiffness matrix/list of nodenumber and coressponding DoFs/dof per nodes
    K = deepcopy(S)
    droptol!(K, 1e-20)
    n, m = size(K) # assume quadratic matrix

    correspondsTogDoF = CoherentStructures.BCTable(ctx, bdata) # new "global" Dofs which are equal to Dof distribution for 1 Dof per node

    new_n = Int(length(unique(correspondsTogDoF))*ndof)
    
    vals = nonzeros(K)
    rows = rowvals(K)

    # Make an empty sparse matrix
    I = Int[]
    sizehint!(I, length(rows))
    J = Int[]
    sizehint!(J, length(rows))
    vals = nonzeros(K)
    V = Float64[]
    sizehint!(V, length(rows))

    convinv = p -> ctx.node_to_dof[p] # converts node number to gDof

    @inbounds for j in 1:m # j: DoF number
        n1, fd1 = nodedoflist[j,:]
        @inbounds for i in nzrange(K,j) # i: DoF number 
            row = rows[i] # welche Reihe
            n2, fd2 = nodedoflist[row,:]
            push!(I, (correspondsTogDoF[convinv(n2)]-1)*ndof+fd2) 
            push!(J, (correspondsTogDoF[convinv(n1)]-1)*ndof+fd1)   
            push!(V, vals[i])
        end
    end
    
    K_red = sparse(I, J, V, new_n, new_n)
    return K_red
end

function applyPBC_r(ctx::CoherentStructures.GridContext, bdata::CoherentStructures.BoundaryData, r::Array{Float64,1}, nodedoflist::Array{Int64,2}, ndof::Int) # gridcontext/boundary data/residual/list of nodenumber and coressponding DoFs/dof per nodes
    n = length(r) 
    correspondsTogDoF = CoherentStructures.BCTable(ctx, bdata) # new "global" Dofs which are equal to Dof distribution for 1 Dof per node
    
    new_n = length(unique(correspondsTogDoF))
    r_red = zeros(Int(new_n*ndof))
    
    convinv = p -> ctx.node_to_dof[p] # converts node number to gDof
    @inbounds for i in 1:n # iterate through all Dofs of old system
        n1, fd1 = nodedoflist[i,:]
        r_red[((correspondsTogDoF[convinv(n1)]-1)*ndof+fd1)] += r[i]
    end
    return r_red
end

function undoPBC(ctx::CoherentStructures.GridContext, bdata::CoherentStructures.BoundaryData, n_dofs::Int, n_dofs_n::Int, ΔΔu_red::Array{Float64,1}, nodedoflist::Array{Int64,2})
            ΔΔu = zeros(n_dofs)
            ΔΔu_r_dofs = [[zeros(nBCDofs(ctx, bdata)) for _ in 1:n_dofs_n]] 
            ΔΔu_dofs = [[zeros(Int(n_dofs/n_dofs_n)) for _ in 1:n_dofs_n]]
            
            # split reduced solver vector into different field dofs
            @inbounds for i in 1:nBCDofs(ctx, bdata) # loop gDof
                @inbounds for j in 1:n_dofs_n
                    ΔΔu_r_dofs[1][j][i] = ΔΔu_red[Int(n_dofs_n*i-(n_dofs_n-j))]   
                end
            end
            
            # undo PBC for each field seperately
            @inbounds for j in 1:n_dofs_n
                ΔΔu_dofs[1][j][:] = undoBCS(ctx, ΔΔu_r_dofs[1][j][:], bdata)
            end
        
            # calculate non reduced solver vector
            convinv = p -> ctx.node_to_dof[p] # converts node number to gDof number 
            @inbounds for k in 1:n_dofs # k: DoF-number
                i, j = nodedoflist[k,:]
                ΔΔu[k] =  ΔΔu_dofs[1][j][convinv(i)]
            end
    return ΔΔu
end