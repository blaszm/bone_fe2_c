include("C:/Program Files (x86)/Gmsh/gmsh-4.5.6-Windows64-sdk/lib/gmsh.jl");

function create_rve_msh(a::Float64, b::Float64, n::Int)
 
    work_dir = pwd()
    cd("mesh")
    
    b_2 = b/2
    
    lc = b/n
    #lcb = b/n^2
    
    gmsh.initialize()

    gmsh.option.setNumber("General.Terminal", 1)

    gmsh.model.add("rve") #new model + name

    # points (0D)
    gmsh.model.geo.addPoint(-a-b_2, -a-b_2, -a-b_2, lc, 1) # x/y/z-Koord, meshSize, tag ("Punktnummer/-name")
    gmsh.model.geo.addPoint(-b_2, -a-b_2,  -a-b_2, lc, 2)
    gmsh.model.geo.addPoint(b_2, -a-b_2, -a-b_2, lc, 3)
    gmsh.model.geo.addPoint(a+b_2, -a-b_2, -a-b_2, lc, 4)
    gmsh.model.geo.addPoint(a+b_2, -b_2, -a-b_2, lc, 5)
    gmsh.model.geo.addPoint(b_2, -b_2, -a-b_2, lc, 6)
    gmsh.model.geo.addPoint(-b_2, -b_2, -a-b_2, lc, 7)
    gmsh.model.geo.addPoint(-a-b_2, -b_2, -a-b_2, lc, 8)
    gmsh.model.geo.addPoint(-a-b_2, b_2, -a-b_2, lc, 9)
    gmsh.model.geo.addPoint(-b_2, b_2, -a-b_2, lc, 10)
    gmsh.model.geo.addPoint(b_2, b_2, -a-b_2, lc, 11)
    gmsh.model.geo.addPoint(a+b_2, b_2, -a-b_2, lc, 12)
    gmsh.model.geo.addPoint(a+b_2, a+b_2, -a-b_2, lc, 13)
    gmsh.model.geo.addPoint(b_2, a+b_2, -a-b_2, lc, 14)
    gmsh.model.geo.addPoint(-b_2, a+b_2, -a-b_2, lc, 15)
    gmsh.model.geo.addPoint(-a-b_2, a+b_2, -a-b_2, lc, 16)

    # lines (1D)
    # horizontal
    gmsh.model.geo.addLine(1, 2, 1) # line between start tag point / en tag point / 3rd argument: line tag
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    
    gmsh.model.geo.addLine(8, 7, 4)
    gmsh.model.geo.addLine(7, 6, 5)
    gmsh.model.geo.addLine(6, 5, 6)
    
    gmsh.model.geo.addLine(9, 10, 7)
    gmsh.model.geo.addLine(10, 11, 8)
    gmsh.model.geo.addLine(11, 12, 9)
    
    gmsh.model.geo.addLine(16, 15, 10)
    gmsh.model.geo.addLine(15, 14, 11)
    gmsh.model.geo.addLine(14, 13, 12)

    # vertikal
    gmsh.model.geo.addLine(1, 8, 13)
    gmsh.model.geo.addLine(8, 9, 14)
    gmsh.model.geo.addLine(9, 16, 15)
    
    gmsh.model.geo.addLine(2, 7, 16)
    gmsh.model.geo.addLine(7, 10, 17)
    gmsh.model.geo.addLine(10, 15, 18)
    
    gmsh.model.geo.addLine(3, 6, 19)
    gmsh.model.geo.addLine(6, 11, 20)
    gmsh.model.geo.addLine(11, 14, 21)
    
    gmsh.model.geo.addLine(4, 5, 22)
    gmsh.model.geo.addLine(5, 12, 23)
    gmsh.model.geo.addLine(12, 13, 24)
    
    
    # surfaces (2D)
    gmsh.model.geo.addCurveLoop([1, 16, -4, -13], 1) # curve tags + Reihenfolge! / loop tag
    gmsh.model.geo.addPlaneSurface([1], 1) # wire (loop) tag / surface tag, ggbf. add. loops for holes

    gmsh.model.geo.addCurveLoop([2, 19, -5, -16], 2)
    gmsh.model.geo.addPlaneSurface([2], 2) 
    
    gmsh.model.geo.addCurveLoop([3, 22, -6, -19], 3)
    gmsh.model.geo.addPlaneSurface([3], 3) 
    
    gmsh.model.geo.addCurveLoop([4, 17, -7, -14], 4)
    gmsh.model.geo.addPlaneSurface([4], 4) 
    
    gmsh.model.geo.addCurveLoop([5, 20, -8, -17], 5)
    gmsh.model.geo.addPlaneSurface([5], 5) 
    
    gmsh.model.geo.addCurveLoop([6, 23, -9, -20], 6)
    gmsh.model.geo.addPlaneSurface([6], 6) 
    
    gmsh.model.geo.addCurveLoop([7, 18, -10, -15], 7)
    gmsh.model.geo.addPlaneSurface([7], 7) 
    
    gmsh.model.geo.addCurveLoop([8, 21, -11, -18], 8)
    gmsh.model.geo.addPlaneSurface([8], 8) 
    
    gmsh.model.geo.addCurveLoop([9, 24, -12, -21], 9)
    gmsh.model.geo.addPlaneSurface([9], 9) 

    # für quad element mesh nötig: Transfinite + Recombine
    for i in 1:9
        gmsh.model.geo.mesh.setTransfiniteSurface(i, "Left") #tag/Ausrichtung
        gmsh.model.geo.mesh.setRecombine(2, i) #dim/tag
    end

    # volumes from extrude (3D) - Layer 1
    DimTags1 = gmsh.model.geo.extrude([(2, 1)], 0, 0, a, [n], [1], true) # dimtags: dim/tag Kombi aller erzeugten Flächen / Volumen, Volumen ist immer an 2. Stelle
    DimTags2 = gmsh.model.geo.extrude([(2, 2)], 0, 0, a, [n], [1], true) # Fläche, die extrudiert werden soll, x/y/z-Länge der Operation, Anzahl Layer, %-
    DimTags3 = gmsh.model.geo.extrude([(2, 3)], 0, 0, a, [n], [1], true)
    DimTags4 = gmsh.model.geo.extrude([(2, 4)], 0, 0, a, [n], [1], true)
    DimTags5 = gmsh.model.geo.extrude([(2, 5)], 0, 0, a, [n], [1], true)
    DimTags6 = gmsh.model.geo.extrude([(2, 6)], 0, 0, a, [n], [1], true)
    DimTags7 = gmsh.model.geo.extrude([(2, 7)], 0, 0, a, [n], [1], true)
    DimTags8 = gmsh.model.geo.extrude([(2, 8)], 0, 0, a, [n], [1], true)
    DimTags9 = gmsh.model.geo.extrude([(2, 9)], 0, 0, a, [n], [1], true)

    # Layer 2
    DimTags10 = gmsh.model.geo.extrude([DimTags1[1]], 0, 0, b, [n], [1], true)
    DimTags11 = gmsh.model.geo.extrude([DimTags2[1]], 0, 0, b, [n], [1], true)
    DimTags12 = gmsh.model.geo.extrude([DimTags3[1]], 0, 0, b, [n], [1], true)
    DimTags13 = gmsh.model.geo.extrude([DimTags4[1]], 0, 0, b, [n], [1], true)
    DimTags14 = gmsh.model.geo.extrude([DimTags5[1]], 0, 0, b, [n], [1], true)
    DimTags15 = gmsh.model.geo.extrude([DimTags6[1]], 0, 0, b, [n], [1], true)
    DimTags16 = gmsh.model.geo.extrude([DimTags7[1]], 0, 0, b, [n], [1], true)
    DimTags17 = gmsh.model.geo.extrude([DimTags8[1]], 0, 0, b, [n], [1], true)
    DimTags18 = gmsh.model.geo.extrude([DimTags9[1]], 0, 0, b, [n], [1], true)
    
    # Layer 3
    DimTags19 = gmsh.model.geo.extrude([DimTags10[1]], 0, 0, a, [n], [1], true)
    DimTags20 = gmsh.model.geo.extrude([DimTags11[1]], 0, 0, a, [n], [1], true)
    DimTags21 = gmsh.model.geo.extrude([DimTags12[1]], 0, 0, a, [n], [1], true)
    DimTags22 = gmsh.model.geo.extrude([DimTags13[1]], 0, 0, a, [n], [1], true)
    DimTags23 = gmsh.model.geo.extrude([DimTags14[1]], 0, 0, a, [n], [1], true)
    DimTags24 = gmsh.model.geo.extrude([DimTags15[1]], 0, 0, a, [n], [1], true)
    DimTags25 = gmsh.model.geo.extrude([DimTags16[1]], 0, 0, a, [n], [1], true)
    DimTags26 = gmsh.model.geo.extrude([DimTags17[1]], 0, 0, a, [n], [1], true)
    DimTags27 = gmsh.model.geo.extrude([DimTags18[1]], 0, 0, a, [n], [1], true)
    
    # Material / physical group
    gmsh.model.addPhysicalGroup(3, [5, 11, 13, 14, 15, 17, 23], 1)
    gmsh.model.setPhysicalName(3, 1, "Material_1")
    gmsh.model.addPhysicalGroup(3, [1, 2, 3, 4, 6, 7, 8, 9, 10, 12, 16, 18, 19, 20, 21, 22, 24, 25, 26, 27], 2)
    gmsh.model.setPhysicalName(3, 2, "Material_2")

    # left right top bottom front back: Flächen sammeln (geht auch in Ferrite direkt sonst)
    

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3) # 3 = dim
    gmsh.write("rve.msh")
    gmsh.finalize()
    
    cd(work_dir)
end # of fct create_rve_msh