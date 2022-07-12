
using SparseArrays

"""
Reads Abaqus-Mesh from *.inp file (supported as export format in GMSH)
    input: file (pasth + name)
    returns: grid, cell_type, [cells,cells_ELSET], [surfaces,surfaces_ELSET]
        grid::Grid -> complete 'Grid' data for Ferrite incl 'boundary_matrix' (with emtpy face- and cell-sets)
        cell_type::AbstractCell -> reader is only suitable for meshes containing one cell type
        cells::Vector{Vector{Integer}} -> global node numbers of each cell
        cells_ELSET::Vector{String} -> Vector with String for each cell containing the corresponding ELSET-name
        surfaces::Vector{Vector{Integer}} -> global node numbers of each surface in order of apperance in imported file
        surfaces_ELSET::Matrix{String} -> Vector with String for each surface face containing the corresponding ELSET-name
"""
function mesh_from_inp(fname::String)

    cell_type = nothing

    # get data from file
    f = open(fname,"r")
    println("  reading mesh file ",fname)
    data = read(f, String)

    # get number of spacial dimensions
    nsd = 0
    volume = findfirst("ELSET=Volume",data)
    surf = findfirst("ELSET=Surface",data)
    if isnothing(volume)
        if isnothing(surf)
            error("Only 2- and 3-dimensional meshes are supported!")
        else
            nsd = 2
        end
    else
        nsd = 3
    end

    # get element type
    ele_type = ""
    if nsd == 3
        posi1 = findlast('=',data[1:volume[1]])
        posi2 = findlast(',',data[1:volume[1]])
    else # nsd ==2
        posi1 = findlast('=',data[1:surf[1]])
        posi2 = findlast(',',data[1:surf[1]])
    end
    ele_type = data[posi1+1:posi2-1]

    # check if element_type is supported and define according surface type
    surf_type = ""
    if nsd == 3
        if ele_type == "C3D4"
            surf_type = "CPS3"
            cell_type = Ferrite.Tetrahedron
        elseif ele_type == "C3D8"
            surf_type = "CPS4"
            cell_type = Ferrite.Hexahedron
        elseif ele_type == "C3D10"
            surf_type = "CPS6"
            cell_type = Ferrite.QuadraticTetrahedron
        else
            error("Non-supported element type ''",ele_type,"'!")
        end
    else # nsd ==2
        if ele_type == "CPS3"
            surf_type = "T3D2"
            cell_type = Ferrite.Triangle
        elseif ele_type == "CPS4"
            surf_type = "T3D2"
            cell_type = Ferrite.Quadrilateral
        elseif ele_type == "CPS6"
            surf_type = "T3D3"
            cell_type = Ferrite.QuadraticTriangle
        else
            error("Non-supported element type ''",ele_type,"'!")
        end
    end
    close(f)

    # get mesh data
    f = open(fname,"r")
    nodes = []
    cells = []
    cells_ELSET = []
    surfaces = []
    surfaces_ELSET = []
    active = ""
    set_name = ""


    while !eof(f) # go line-wise through whole file
        # find header
        line = strip(readline(f))
        if line[1] == '*' # header lines start with '*'
            if occursin("*NODE",line)
                active = "node"
            elseif occursin(ele_type,line)
                active = "cell"
                set_name = line[findlast('=',line)+1:end]
            elseif occursin(surf_type,line)
                active = "surface"
                set_name = line[findlast('=',line)+1:end]
            else
                active = ""
            end

        else # data_lines start with '*'
            line_data = strip.(split(line,',',keepempty=false))
            if active == "node"
                push!(nodes,parse.(Float64, line_data)[2:nsd+1])
            elseif active == "cell"
                push!(cells,parse.(Int, line_data)[2:end])
                push!(cells_ELSET,set_name)
            elseif active == "surface"
                push!(surfaces,parse.(Int, line_data)[2:end])
                push!(surfaces_ELSET,set_name)
            end
        end
    end
    close(f)

    # find hanging nodes
    nnodes = size(nodes,1)
    ncells = size(cells,1)
    adj_cells = get_con_cells(nnodes,ncells,cells)
    hanging_nodes = Integer[]
    for (node,cells) in enumerate(adj_cells)
        if isempty(cells)
            push!(hanging_nodes,node)
        end
    end
    sort!(hanging_nodes)

    # removing hanging nodes from data
    if isempty(hanging_nodes)
        println("  no hanging nodes detected.")
    else
        println("  deleting hanging nodes ",hanging_nodes)
        # delete coordinates
        deleteat!(nodes,hanging_nodes)
        # reduce node numbers in cells and surfaces
        for h_node in reverse(hanging_nodes)
            for cell_nodes in cells
                hits = findall(x -> x > h_node,cell_nodes)
                cell_nodes[hits] = cell_nodes[hits] .- 1
            end
            for surf_nodes in surfaces
                hits = findall(x -> x > h_node,surf_nodes)
                surf_nodes[hits] = surf_nodes[hits] .- 1
            end
        end
    end

    # convert to Ferrite and construct set_grid
    nnodes = size(nodes,1)
    ncells = size(cells,1)
    println("  found-$ele_type mesh with $ncells cells and $nnodes nodes...")
    j_nodes = Vector{Ferrite.Node{nsd,Float64}}(undef,nnodes)
    j_cells = Vector{cell_type}(undef,ncells)
    for i = 1:nnodes
        j_nodes[i] = Ferrite.Node(Vec{nsd}(nodes[i]))
    end
    for i = 1:ncells
        j_cells[i] = cell_type(Tuple(cells[i]))
    end

    # number of faces per cell
    nfc = 0
    if cell_type == Quadrilateral
        nfc = 4
    elseif cell_type == Triangle || cell_type == QuadraticTriangle
        nfc = 3
    elseif cell_type == Hexahedron
        nfc = 6
    elseif cell_type == Tetrahedron || cell_type == QuadraticTetrahedron
        nfc = 4
    end

    # faces on boundary
    boun_mat = Array{Bool,2}(undef, nfc, ncells)
    fill!(boun_mat,false)

    adj_cells = get_con_cells(nnodes,ncells,cells)

    for (i,surf) in enumerate(surfaces)
        found = false
        progress = Integer(round(100.0*i/size(surfaces,1),digits=0))
        #print("\r  converting data for Ferrite... $progress %")

        for cell_id in adj_cells[surf[1]] # go through all cells connected to first node of surf
            cell_faces = Ferrite.faces(j_cells[cell_id])
            for (face_id,face) in enumerate(cell_faces)
                if face ⊆ surf
                    boun_mat[face_id,cell_id] = true
                    found = true
                end
                found ? break : nothing
            end
            found ? break : nothing
        end
    end
    print("\n")
    grid = Ferrite.Grid(j_cells, j_nodes, boundary_matrix=sparse(boun_mat))

    return grid, cell_type, [cells,cells_ELSET], [surfaces,surfaces_ELSET]

end

function get_con_cells(nn::Integer,ne::Integer,cells)
    con_cells = Vector{Vector{Integer}}(undef,nn)
    for i=1:size(con_cells,1)
        con_cells[i] = []
    end
    # filling with data
    for ele = 1:ne # Loop through all element    # Loop through all element
        lNod = cells[ele] # local nodes
        for i = 1:size(lNod,1)                  # loop through all local nodes
            push!(con_cells[lNod[i]],ele)
        end
    end
    return con_cells
end
