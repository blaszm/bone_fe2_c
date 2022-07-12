# Set up cell- and facevalues -> Interpolationsregeln wählen (Grad Ansatzfunktionen für Geometrie und Potentialfelder)
function create_values(etype::DataType) # Aufruf mit linear, d.h. Lagrange{2,RefCube,1}() = interpolation_u
    if(etype == Hexahedron)
        # quadrature rules (Gaußintegration Regel)
        # QuadratureRule{dim,Shape}(Order) : "Normal" 1/sqrt3 ist Order 2
        qr      = QuadratureRule{3,RefCube}(2)
        face_qr = QuadratureRule{2,RefCube}(2)

        # geometric interpolation
        # Lagrange{}(dim,Shape,Order) : "Normal" Ordnung 1, d.h. lineare Ansatzfunktionen 
        interpolation_geom = Lagrange{3,RefCube,1}()

        # cell and facevalues for u
        cellvalues_u = CellVectorValues(qr, Lagrange{3,RefCube,1}(), interpolation_geom)
        cellvalues_phi = CellScalarValues(qr, Lagrange{3,RefCube,1}(), interpolation_geom)
        cellvalues_A = CellVectorValues(qr, Lagrange{3,RefCube,1}(), interpolation_geom)

        facevalues_u = FaceVectorValues(face_qr, Lagrange{3,RefCube,1}(), interpolation_geom) 
    elseif(etype == Tetrahedron)
        qr      = QuadratureRule{3,RefTetrahedron}(2)
        face_qr = QuadratureRule{2,RefTetrahedron}(2)

        # geometric interpolation
        # Lagrange{}(dim,Shape,Order) : "Normal" Ordnung 1, d.h. lineare Ansatzfunktionen 
        interpolation_geom = Lagrange{3,RefTetrahedron,1}()

        # cell and facevalues for u
        cellvalues_u = CellVectorValues(qr, Lagrange{3,RefTetrahedron,1}(), interpolation_geom)
        cellvalues_phi = CellScalarValues(qr, Lagrange{3,RefTetrahedron,1}(), interpolation_geom)
        cellvalues_A = CellVectorValues(qr, Lagrange{3,RefTetrahedron,1}(), interpolation_geom)

        facevalues_u = FaceVectorValues(face_qr, Lagrange{3,RefTetrahedron,1}(), interpolation_geom)
    else
        println("Error: element type not supported for cellvalues")
    end
        
    return cellvalues_u, cellvalues_phi, cellvalues_A, facevalues_u
end;

# DofHandler (welche Freiheitsgrade pro Knoten und welche Namen)
function create_dofhandler(grid::Grid, etype::DataType)
    dh = DofHandler(grid)
    if(etype == Hexahedron)
        push!(dh, :u, 3, Lagrange{3,RefCube,1}()) # mechanical displacement 3 DoF
        push!(dh, :phi, 1, Lagrange{3,RefCube,1}()) # electric scalar potential 1 DoF
        push!(dh, :A, 3, Lagrange{3,RefCube,1}()) # magnetic vector potential 3 DoF
    elseif(etype == Tetrahedron)
        push!(dh, :u, 3, Lagrange{3,RefTetrahedron,1}()) # mechanical displacement 3 DoF
        push!(dh, :phi, 1, Lagrange{3,RefTetrahedron,1}()) # electric scalar potential 1 DoF
        push!(dh, :A, 3, Lagrange{3,RefTetrahedron,1}()) # magnetic vector potential 3 DoF
    else
        println("Error: element type not supported for dofhandler")
    end
    close!(dh)
    return dh
end;

function create_dofhandler_helper(grid::Grid, etype::DataType)
    dh = DofHandler(grid)
    if(etype == Hexahedron)
        push!(dh, :T, 1, Lagrange{3,RefCube,1}())
    elseif(etype == Tetrahedron)
        push!(dh, :T, 1, Lagrange{3,RefTetrahedron,1}())
    else
        println("Error: element type not supported for dofhandler helper")
    end
    close!(dh)
    return dh
end; 

# boundary conditions
function create_bc(dh::DofHandler) 
    dbc = ConstraintHandler(dh)
    # fix corner nodes completely
    add!(dbc, Dirichlet(:u, getnodeset(dh.grid, "corners"), (x,t) -> zero(Vec{3}), [1,2,3]))
    add!(dbc, Dirichlet(:phi, getnodeset(dh.grid, "corners"), (x,t) -> zero(Vec{1}), [1])) 
    add!(dbc, Dirichlet(:A, getnodeset(dh.grid, "corners"), (x,t) -> zero(Vec{3}), [1,2,3])) 
    close!(dbc)
    t = 0.0 # t = time
    Ferrite.update!(dbc, t)
    return dbc
end;