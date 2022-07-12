using Ferrite
using Printf

const FILE_FORMAT = ("4.1", "0", "8") # 2.2 alt 

function verify_next_line(s, f)
    line = readline(f)
    @assert strip(line) == s
end

"""
mesh = read_msh_file("filename.msh")
Read a .msh file created by Gmsh and return the corresponding 
Mesh object.
"""
function read_msh_file(fname::String)
    
    nodenumbers = Int[] # speichert die Nummern der Knoten
    coord = Float64[] # speichert die Koordinaten der Knoten
    nonodes = 0 # Gesamtanzahl der Knoten
    line = ""
    
    elms = String[] # Liste der Elemente und Knotenzugehörigkeiten als string
    elmcount = Array{Int,1}[] # speichert die Elementart und Anzahl der Elemente jedes Blocks
    
    physdim  = Dict{String,Int}() # dictionaries für das Sammeln der physischen Bereiche
    physnum  = Dict{String,Int}()
    physname = Dict{Int,String}()
    
    physentities = Array{Int,1}[] # dim entity / entity no / phys. name no of entity
    elmentities = Array{Int,1}[] # element no / entity no

    f = open(fname,"r") 
    while !eof(f) # schaut, ob Ende des Files erreicht wurde
        
        header = strip( readline(f) ) # prüft, in welchem Bereich wir sind
        if header == "\$MeshFormat" # prüft, ob MeshFormat zur aktuellen Vers passt (habe ich manuell geändert, dadurch sind Fehler entstanden)
            fmt = split( readline(f) )
            @inbounds for j = 1:3
                @assert fmt[j] == FILE_FORMAT[j]
            end
            verify_next_line("\$EndMeshFormat", f)
            
        elseif header == "\$PhysicalNames" 
            line = readline(f)
            nonames = parse(Int, line)
            @inbounds for n = 1:nonames
                line = readline(f)
                s = split(line)
                dimen = parse(Int, s[1])
                num   = parse(Int, s[2])
                name  = strip(s[3], '"')
                physdim[name] = dimen
                physnum[name] = num
                physname[num] = name
            end
            verify_next_line("\$EndPhysicalNames", f)       
         
        
        elseif header == "\$Entities" # hier ist Materialzugehörigkeit versteckt! (Stelle 9 der 2d/3d Zeilen)
            line = readline(f)
            s = split(line)
            noent_0d = parse(Int, s[1]) 
            noent_1d = parse(Int, s[2])
            noent_2d = parse(Int, s[3])
            noent_3d = parse(Int, s[4])
            exist_2d = 2 in values(physdim) ? true : false
            exist_3d = 3 in values(physdim) ? true : false
            @inbounds for i in 1:(noent_0d+noent_1d)
                readline(f)
            end
            @inbounds for i in 1:noent_2d
                line = readline(f)
                if(exist_2d)
                    s = split(line)
                    a = zeros(3)
                    a[1] = 2 # dim
                    a[2] = parse(Int, s[1]) # entity
                    a[3] = parse(Int, s[9]) # phys name no
                    push!(physentities, a)
                end
            end
            @inbounds for i in 1:noent_3d
                line = readline(f)
                if(exist_3d)
                    s = split(line)
                    a = zeros(3)
                    a[1] = 3 # dim 
                    a[2] = parse(Int, s[1]) # entity
                    a[3] = parse(Int, s[9]) # phys name no
                    push!(physentities, a)
                end
            end
            verify_next_line("\$EndEntities", f)  
        
            
        elseif header == "\$Nodes"
            line = readline(f)
            s = split(line)
            nodeblocks = parse(Int, s[1])
            nonodes = parse(Int, s[2])
            @inbounds for b = 1:nodeblocks # Anzahl Blöcke
                line = readline(f)
                s = split(line)
                nonodescurrblock = parse(Int, s[4]) # Anzahl Nodes im aktuellen Block
                @inbounds for n = 1:nonodescurrblock
                    line = readline(f)
                    push!(nodenumbers, parse(Int, line)) # Knotennummern abspeichern
                end
                @inbounds for n = 1:nonodescurrblock
                    line = readline(f)
                    s = split(line)
                    @inbounds for j = 1:3
                       push!(coord, round(parse(Float64, s[j]), digits=8)) 
                    end
                end    
            end
            verify_next_line("\$EndNodes", f) 
            
        elseif header == "\$Elements" 
            line = readline(f)
            s = split(line)
            elemblocks = parse(Int, s[1])
            noelems = parse(Int, s[2])
            el = 1
            @inbounds for b in 1:elemblocks
                line = readline(f)
                s = split(line)
                entityno = parse(Int, s[2]) #!! integer 1,..,max
                elmtype = parse(Int, s[3])
                noelemscurrblock = parse(Int, s[4])
                push!(elmcount, [elmtype, noelemscurrblock])
                @inbounds for n = 1:noelemscurrblock
                    line = readline(f)
                    push!(elms, line)
                    a = zeros(2)
                    a[1] = el
                    a[2] = entityno
                    push!(elmentities, a)
                    el += 1
                end
            end
            verify_next_line("\$EndElements", f)
            
        elseif header == "\$Comment" 
            while true
                line = readline(f)
                if strip(line) == "\$EndComment"
                    break
                end
            end
        else
            msg = string("Unknown section header: ", line)
            error(msg)
        end
    end
    close(f)
    
# ---- ab hier Umwandlung eingelesenes Mesh -> Ferrite.mesh
    
    # element type
    # vorerst nur ein Elementtyp zugelassen
    oneelmtype = true
    @inbounds for i in 2:length(elmcount)
        if(elmcount[1][1] != elmcount[i][1])
            oneelmtype = false
            break
        end
    end
    
    @assert oneelmtype
    
    gmsh_elmtype_key = elmcount[1][1]
    
    # nodes
    coord = reshape(coord, (3,nonodes)) # Vorsicht hier alternative für 2d einbauen!
    
    @assert length(nodenumbers) == nonodes
    @assert issorted(nodenumbers) # Achtung: Das klappt so nur, wenn die Knotennummern schon sortiert sind
    # Wäre das nicht der Fall: Nicht nur die Knotennummern, sondern auch die Koordinaten müssten sortiert werden! (ggnf. mit dictionary oder so möglich??)
    
    if(gmsh_elmtype_key == 3) # 2d-Quadrilateral elements
        nodes = Vector{Ferrite.Node{2,Float64}}(undef, nonodes)
        @inbounds for i in 1:nonodes
            nodes[i] = Ferrite.Node(Vec{2}(coord[1:2,i]))
        end
    elseif(gmsh_elmtype_key == 4) # 3d-Tetrahedron elements
        nodes = Vector{Ferrite.Node{3,Float64}}(undef, nonodes)
        @inbounds for i in 1:nonodes
            nodes[i] = Ferrite.Node(Vec{3}(coord[:,i]))
        end
    elseif(gmsh_elmtype_key == 5) # 3d-Hexahedron elements
        nodes = Vector{Ferrite.Node{3,Float64}}(undef, nonodes)
        @inbounds for i in 1:nonodes
            nodes[i] = Ferrite.Node(Vec{3}(coord[:,i]))
        end
    end
        
    # elements (altenativ mit cases)
    #if(gmsh_elmtype_key==1) # line
    #    cells = Vector{Ferrite.Line}(undef, length(elms))
    #elseif(gmsh_elmtype_key==2)
    #    cells = Vector{Ferrite.Triangle}(undef, length(elms))
    #elseif(gmsh_elmtype_key==3)
    if(gmsh_elmtype_key==3)
        cells = Vector{Ferrite.Quadrilateral}(undef, length(elms))
        @inbounds for (i, elm) in enumerate(elms)
            s = split(elm)
            cellnodes = (parse(Int, s[2]), parse(Int, s[3]), parse(Int, s[4]), parse(Int, s[5]))
            cells[i] = Ferrite.Quadrilateral(cellnodes)
        end
    elseif(gmsh_elmtype_key==4)
        cells = Vector{Ferrite.Tetrahedron}(undef, length(elms))
        @inbounds for (i, elm) in enumerate(elms)
            s = split(elm)
            cellnodes = (parse(Int, s[2]), parse(Int, s[3]), parse(Int, s[4]), parse(Int, s[5]))
            cells[i] = Ferrite.Tetrahedron(cellnodes)
        end
    elseif(gmsh_elmtype_key==5)
        cells = Vector{Ferrite.Hexahedron}(undef, length(elms))
        @inbounds for (i, elm) in enumerate(elms)
            s = split(elm)
            cellnodes = (parse(Int, s[2]), parse(Int, s[3]), parse(Int, s[4]), parse(Int, s[5]), parse(Int, s[6]), parse(Int, s[7]), parse(Int, s[8]), parse(Int, s[9]))
            cells[i] = Ferrite.Hexahedron(cellnodes)
        end
    end

#----------------------------
    
    # cellsets
    cellsets = Dict{String, Set{Int}}()
    
    if(gmsh_elmtype_key==3) #2d
        @inbounds for i in 1:length(physname) # loop physical names
            name = physname[i]
            if(physdim[name] == 2)
                number = physnum[name]
                cellsofset = Set{Int}()
                @inbounds for j in 1:length(physentities) # loop entities
                    if(physentities[j][3]==number) # falls das erfüllt ist: finde alle Elemente aus dieser entity und füge sie hinzu
                        entityno = physentities[j][2]
                        @inbounds for k in 1:length(elmentities)
                            if(elmentities[k][2] == entityno)
                                new_cell = Set{Int}(elmentities[k][1])
                                cellsofset = union(cellsofset, new_cell)
                            end # if elementities
                        end # for elmentities
                    end
                end # for physentities
                cellsets[name] = cellsofset
            end
        end
    elseif(gmsh_elmtype_key==4 || gmsh_elmtype_key==5) #3d
        @inbounds for i in 1:length(physname) # loop physical names
            name = physname[i]
            if(physdim[name] == 3)
                number = physnum[name]
                cellsofset = Set{Int}()
                @inbounds for j in 1:length(physentities) # loop entities
                    if(physentities[j][3]==number) # falls das erfüllt ist: finde alle Elemente aus dieser entity und füge sie hinzu
                        entityno = physentities[j][2]
                        @inbounds for k in 1:length(elmentities)
                            if(elmentities[k][2] == entityno) 
                                new_cell = Set{Int}(elmentities[k][1])
                                cellsofset = union(cellsofset, new_cell)
                            end # if elementities
                        end # for elmentities
                    end
                end # for physentities
                cellsets[name] = cellsofset
            end
        end
    end
    
    return Ferrite.Grid(cells, nodes; cellsets = cellsets) 
end