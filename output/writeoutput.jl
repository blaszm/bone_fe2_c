function WriteVTK.vtk_grid(filename::AbstractString, grid::Grid{dim,C,T}; compress=false) where {dim,C,T} # custom overload with kwargs
    cls = MeshCell[]
    for cell in grid.cells
        celltype = Ferrite.cell_to_vtkcell(typeof(cell))
        push!(cls, MeshCell(celltype, collect(cell.nodes)))
    end
    coords = reshape(reinterpret(T, getnodes(grid)), (dim, getnnodes(grid)))
	if(compress)
		return vtk_grid(filename, coords, cls; compress=true)
	else
		return vtk_grid(filename, coords, cls; compress=false) # MB: added compress option
	end
end

function writeoutput(cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, dh::DofHandler, dbc::ConstraintHandler, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}, grid::Grid, t::Int, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, material_number::Array{Float64,1}, globgpnumber::Int, pvd::WriteVTK.CollectionFile, etype_micro::DataType) where {dim} 
    
    work_dir = pwd()
    cd("Results")
    
    # Output stresses / strain / inelastic strain - x,y,z,xy,yz,xz
    sigma_xx = collect_stresses(cellvalues_u, dh, states, 1)
    sigma_yy = collect_stresses(cellvalues_u, dh, states, 2)
    sigma_zz = collect_stresses(cellvalues_u, dh, states, 3)
    sigma_xy = collect_stresses(cellvalues_u, dh, states, 4)
    sigma_yz = collect_stresses(cellvalues_u, dh, states, 5)
    sigma_xz = collect_stresses(cellvalues_u, dh, states, 6)
    eps_xx = collect_strain(cellvalues_u, dh, states, 1)
    eps_yy = collect_strain(cellvalues_u, dh, states, 2)
    eps_zz = collect_strain(cellvalues_u, dh, states, 3)
    eps_xy = collect_strain(cellvalues_u, dh, states, 4)
    eps_yz = collect_strain(cellvalues_u, dh, states, 5)
    eps_xz = collect_strain(cellvalues_u, dh, states, 6)
    epsi_xx = collect_inelstrain(cellvalues_u, dh, states, 1)
    epsi_yy = collect_inelstrain(cellvalues_u, dh, states, 2)
    epsi_zz = collect_inelstrain(cellvalues_u, dh, states, 3)
    epsi_xy = collect_inelstrain(cellvalues_u, dh, states, 4)
    epsi_yz = collect_inelstrain(cellvalues_u, dh, states, 5)
    epsi_xz = collect_inelstrain(cellvalues_u, dh, states, 6)
       
    # Output other fields
    D_out = collect_D(cellvalues_u, dh, states)
    Dp_out = collect_Dp(cellvalues_u, dh, states)
    E_out = collect_E(cellvalues_u, dh, states)
    J_out = collect_J(cellvalues_u, dh, states)
    H_out = collect_H(cellvalues_u, dh, states)
    B_out = collect_B(cellvalues_u, dh, states)
        
    if(etype_micro == Hexahedron)
        projector = L2Projector(cellvalues_phi, Lagrange{3,RefCube,1}(), grid)
    elseif(etype_micro == Tetrahedron)
        projector = L2Projector(cellvalues_phi, Lagrange{3,RefTetrahedron,1}(), grid)
    else
        println("Error: unsupported element type for output")
    end
     
    sigma_xx_nodes = project(sigma_xx, projector)
    sigma_yy_nodes = project(sigma_yy, projector)
    sigma_zz_nodes = project(sigma_zz, projector)
    sigma_xy_nodes = project(sigma_xy, projector)
    sigma_yz_nodes = project(sigma_yz, projector)
    sigma_xz_nodes = project(sigma_xz, projector)
    eps_xx_nodes = project(eps_xx, projector)
    eps_yy_nodes = project(eps_yy, projector)
    eps_zz_nodes = project(eps_zz, projector)
    eps_xy_nodes = project(eps_xy, projector)
    eps_yz_nodes = project(eps_yz, projector)
    eps_xz_nodes = project(eps_xz, projector)
    epsi_xx_nodes = project(epsi_xx, projector)
    epsi_yy_nodes = project(epsi_yy, projector)
    epsi_zz_nodes = project(epsi_zz, projector)
    epsi_xy_nodes = project(epsi_xy, projector)
    epsi_yz_nodes = project(epsi_yz, projector)
    epsi_xz_nodes = project(epsi_xz, projector)
        
    D_nodes = project(D_out, projector)
    Dp_nodes = project(Dp_out, projector)
    E_nodes = project(E_out, projector)
    J_nodes = project(J_out, projector)
    H_nodes = project(H_out, projector)
    B_nodes = project(B_out, projector)
        
    time = @sprintf("%3.3i", t)
    gpnumber = @sprintf("%5.5i", globgpnumber)
    filename = string("micro_", gpnumber, "_", time)
       
    vtk_grid(filename, grid; compress=false) do vtk # erstellt Paraview-Output
        vtk_point_data(vtk, dh, u) 
        vtk_point_data(vtk, dh, v, "p") 
        vtk_point_data(vtk, dh, a, "pp") 
        
        vtk_point_data(vtk, sigma_xx_nodes, "sigma_xx")
        vtk_point_data(vtk, sigma_yy_nodes, "sigma_yy")
        vtk_point_data(vtk, sigma_zz_nodes, "sigma_zz")
        vtk_point_data(vtk, sigma_xy_nodes, "sigma_xy")
        vtk_point_data(vtk, sigma_yz_nodes, "sigma_yz")
        vtk_point_data(vtk, sigma_xz_nodes, "sigma_xz")
            
        vtk_point_data(vtk, eps_xx_nodes, "eps_xx")
        vtk_point_data(vtk, eps_yy_nodes, "eps_yy")
        vtk_point_data(vtk, eps_zz_nodes, "eps_zz")
        vtk_point_data(vtk, eps_xy_nodes, "eps_xy")
        vtk_point_data(vtk, eps_yz_nodes, "eps_yz")
        vtk_point_data(vtk, eps_xz_nodes, "eps_xz")
            
        vtk_point_data(vtk, epsi_xx_nodes, "epsi_xx")
        vtk_point_data(vtk, epsi_yy_nodes, "epsi_yy")
        vtk_point_data(vtk, epsi_zz_nodes, "epsi_zz")
        vtk_point_data(vtk, epsi_xy_nodes, "epsi_xy")
        vtk_point_data(vtk, epsi_yz_nodes, "epsi_yz")
        vtk_point_data(vtk, epsi_xz_nodes, "epsi_xz")
            
        vtk_point_data(vtk, D_nodes, "D")
        vtk_point_data(vtk, Dp_nodes, "Dp")
        vtk_point_data(vtk, E_nodes, "E")
        vtk_point_data(vtk, J_nodes, "J")
        vtk_point_data(vtk, H_nodes, "H")
        vtk_point_data(vtk, B_nodes, "B")
            
        vtk_cell_data(vtk, material_number, "Mat")
        pvd[t-1] = vtk
    end # of export
    cd(work_dir)
    return
end

function writebc(dbc, grid::Grid, filename::String)
    work_dir = pwd()
    cd("Results")
    vtk_grid(filename, grid; compress=false) do vtk
        vtk_point_data(vtk, dbc)
    end
    cd(work_dir)
    return
end

function initpvd(globgpnumber)
    work_dir = pwd()
    cd("Results")
    gpnumber = @sprintf("%5.5i", globgpnumber)
    filename = string("pvd_micro_", gpnumber)
    pvd = paraview_collection(filename)
    cd(work_dir)
    return pvd
end

function savepvd(pvd)
    work_dir = pwd()
    cd("Results")
    vtk_save(pvd)
    cd(work_dir)
    return
end

function writeoutput_macro(cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, dh::DofHandler, dbc::ConstraintHandler, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}, grid::Grid, t::Int, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, pvd::WriteVTK.CollectionFile, r::Array{Float64,1}, etype_macro::DataType) where {dim}
    
    work_dir = pwd()
    cd("Results")
    
    # Output stresses / strain / inelastic strain - x,y,z,xy,yz,xz
    sigma_xx = collect_stresses(cellvalues_u, dh, states, 1)
    sigma_yy = collect_stresses(cellvalues_u, dh, states, 2)
    sigma_zz = collect_stresses(cellvalues_u, dh, states, 3)
    sigma_xy = collect_stresses(cellvalues_u, dh, states, 4)
    sigma_yz = collect_stresses(cellvalues_u, dh, states, 5)
    sigma_xz = collect_stresses(cellvalues_u, dh, states, 6)
    eps_xx = collect_strain(cellvalues_u, dh, states, 1)
    eps_yy = collect_strain(cellvalues_u, dh, states, 2)
    eps_zz = collect_strain(cellvalues_u, dh, states, 3)
    eps_xy = collect_strain(cellvalues_u, dh, states, 4)
    eps_yz = collect_strain(cellvalues_u, dh, states, 5)
    eps_xz = collect_strain(cellvalues_u, dh, states, 6)
       
    # Output other fields
    D_out = collect_D(cellvalues_u, dh, states)
    Dp_out = collect_Dp(cellvalues_u, dh, states)
    E_out = collect_E(cellvalues_u, dh, states)
    J_out = collect_J(cellvalues_u, dh, states)
    H_out = collect_H(cellvalues_u, dh, states)
    B_out = collect_B(cellvalues_u, dh, states)
        
    if(etype_macro == Hexahedron)
        projector = L2Projector(cellvalues_phi, Lagrange{3,RefCube,1}(), grid)
    elseif(etype_macro == Tetrahedron)
        projector = L2Projector(cellvalues_phi, Lagrange{3,RefTetrahedron,1}(), grid)
    else
        println("Error: unsupported element type for output")
    end
    
    sigma_xx_nodes = project(sigma_xx, projector)
    sigma_yy_nodes = project(sigma_yy, projector)
    sigma_zz_nodes = project(sigma_zz, projector)
    sigma_xy_nodes = project(sigma_xy, projector)
    sigma_yz_nodes = project(sigma_yz, projector)
    sigma_xz_nodes = project(sigma_xz, projector)
    eps_xx_nodes = project(eps_xx, projector)
    eps_yy_nodes = project(eps_yy, projector)
    eps_zz_nodes = project(eps_zz, projector)
    eps_xy_nodes = project(eps_xy, projector)
    eps_yz_nodes = project(eps_yz, projector)
    eps_xz_nodes = project(eps_xz, projector)
        
    D_nodes = project(D_out, projector)
    Dp_nodes = project(Dp_out, projector)
    E_nodes = project(E_out, projector)
    J_nodes = project(J_out, projector)
    H_nodes = project(H_out, projector)
    B_nodes = project(B_out, projector)
        
    time = @sprintf("%3.3i", t)
    filename = string("macro_", time)
       
    vtk_grid(filename, grid; compress=false) do vtk # erstellt Paraview-Output
        vtk_point_data(vtk, dh, u) 
        vtk_point_data(vtk, dh, v, "p") 
        vtk_point_data(vtk, dh, a, "pp") 
        
        vtk_point_data(vtk, sigma_xx_nodes, "sigma_xx")
        vtk_point_data(vtk, sigma_yy_nodes, "sigma_yy")
        vtk_point_data(vtk, sigma_zz_nodes, "sigma_zz")
        vtk_point_data(vtk, sigma_xy_nodes, "sigma_xy")
        vtk_point_data(vtk, sigma_yz_nodes, "sigma_yz")
        vtk_point_data(vtk, sigma_xz_nodes, "sigma_xz")
            
        vtk_point_data(vtk, eps_xx_nodes, "eps_xx")
        vtk_point_data(vtk, eps_yy_nodes, "eps_yy")
        vtk_point_data(vtk, eps_zz_nodes, "eps_zz")
        vtk_point_data(vtk, eps_xy_nodes, "eps_xy")
        vtk_point_data(vtk, eps_yz_nodes, "eps_yz")
        vtk_point_data(vtk, eps_xz_nodes, "eps_xz")

        vtk_point_data(vtk, D_nodes, "D")
        vtk_point_data(vtk, Dp_nodes, "Dp")
        vtk_point_data(vtk, E_nodes, "E")
        vtk_point_data(vtk, J_nodes, "J")
        vtk_point_data(vtk, H_nodes, "H")
        vtk_point_data(vtk, B_nodes, "B")
            
        pvd[t-1] = vtk
    end # of export
    
    cd(work_dir)
    return
end

function initpvd_macro()
    work_dir = pwd()
    cd("Results")
    filename = "pvd_macro"
    pvd = paraview_collection(filename)
    cd(work_dir)
    return pvd
end

