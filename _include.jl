using Ferrite
using CoherentStructures, Distances
using BlockArrays, SparseArrays, LinearAlgebra, Tensors
using TimerOutputs, Printf, HDF5
using KrylovMethods
using Distributed

#include("mesh/rvemesh_constr.jl") # mesh creator
include("mesh/geometry_check.jl") # Dustin BC helper
include("mesh/inp_reader.jl") # Dustin inp/Abaqus reader
include("mesh/gmsh_reader.jl") # converts .msh file to Ferrite grid object
include("mesh/mesh.jl"); # mesh tools micro system
include("macroscale/material_macro.jl") # material parameters macro system
include("material/material.jl"); # material modell micro system
include("macroscale/tangent_macro.jl")
include("material/assemble.jl") # assembly micro system
include("material/elmt01.jl") # element routine cortical bone
include("material/elmt02.jl") # element routine bone marrow
include("material/voigt.jl") # tools to convert quantities to Voigt notation
include("material/volaver.jl") # calculate volume averages of quantities
include("material/states.jl") # output state and flux quantities from states
include("material/applypbc.jl") # periodic boundary conditions helper
include("material/solve_RVE.jl") # micro modell solver
include("macroscale/mesh_macro.jl") # mesh tools macro system
include("macroscale/assemble_macro.jl") # assembly macro system
include("macroscale/elmt_macro.jl") # element routine macro system
include("macroscale/model_macro.jl") # different macro models
include("output/writeinfo.jl") # output info
include("output/writeoutput.jl") # output Paraview
include("output/hdffiles.jl") # read, write and copy hdf-files