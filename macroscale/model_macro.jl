function load_model(modeltype::String; option::String="default")
    path = string("mesh/", modeltype,".inp")
    grid, cell_type, cellset, faceset = mesh_from_inp(path)
    
    if(modeltype == "cylinder")
        if(option == "default")
            addfaceset!(grid, "leftsurf", x -> x[3] ≈ 0.0)
            addfaceset!(grid, "rightsurf", x -> x[3] ≈ 0.3)
            addfaceset!(grid, "middlesurf", difference(in_rectangle([-0.015,0.015,0.15],[0.04,0.04,0.05]),in_cylinder([0.0,0.0,0.15],[0.0,0.0,1.0],0.013,0.3)))
            addfaceset!(grid, "leftgrounding", x -> norm(x[1]) < 0.005 && norm(x[2]) < 0.005 && x[3] ≈ 0.0)
        elseif(option == "torsion")
            addfaceset!(grid, "leftsurf", x -> x[3] ≈ 0.0)       
            addfaceset!(grid, "leftgrounding", x -> norm(x[1]) < 0.005 && norm(x[2]) < 0.005 && x[3] ≈ 0.0)
            addnodeset!(grid, "torstop", x -> x[1] ≈ 0.015 && x[2] ≈ 0.0 && x[3] ≈ 0.3)
            addnodeset!(grid, "torsbottom", x -> x[1] ≈ -0.015 && x[2] ≈ 0.0 && x[3] ≈ 0.3)
        end
    elseif(modeltype == "femurbonescaled")
        #addnodeset!(grid, "leftsinglenode", x -> x[1] ≈ -0.102504 && x[2] ≈ -0.0645889 && x[3] ≈ 0.408411)

        addfaceset!(grid, "leftsinglenode", in_sphere([-0.08,-0.06459,0.4084], 0.01))
        leftsurface = getsurfaceset(grid,grid.facesets["leftsinglenode"])
        addfaceset!(grid, "leftgrounding", leftsurface)
        
        addfaceset!(grid, "left", in_sphere([-0.07,-0.055,0.39], 0.045))
        leftsurface = getsurfaceset(grid,grid.facesets["left"])
        addfaceset!(grid, "leftsurf", leftsurface)

        addfaceset!(grid, "right", in_sphere([-0.105,-0.078,0.84], 0.045))
        rightsurface = getsurfaceset(grid,grid.facesets["right"])
        addfaceset!(grid, "rightsurf", rightsurface)

        addfaceset!(grid, "middle", in_sphere([-0.06,-0.09,0.61], 0.045))
        middlesurface = getsurfaceset(grid,grid.facesets["middle"])
        addfaceset!(grid, "middlesurf", middlesurface)
    else
        println("Error: macro model not found")
    end
    
    return grid, cell_type, cellset, faceset
end