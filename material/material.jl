# --- cortical bone
struct CorticalBone{T, S <: AbstractArray{T, 2}} # T = real (number type)
    E::T # Youngs modulus
    ν::T # Poisson ratio
    Cᵉ::S # Elastic stiffness tensor
    ϵ₁::T # dielectric constant 11
    ϵ₃::T # dielectric constant 33
    ϵₜ::S #  dielectric tensor
    μᶜ::T # inverse of magnetic constant
    μⁱ::S # inverse of magnetic tensor 
    eₚ14::T # dielectric constant 14
    eₚ::S # piezoelectric tensor
end

function CorticalBone(E::Float64, ν::Float64, ϵ₁::Float64, ϵ₃::Float64, μᶜ::Float64, eₚ14::Float64) # Constructor zum struct
    # elastic tensor
    Cᵉ = zeros(6,6)
    
    Cᵉ[1,1] = 1 - ν
    Cᵉ[1,2] = ν
    Cᵉ[1,3] = ν
    
    Cᵉ[2,1] = ν
    Cᵉ[2,2] = 1 - ν
    Cᵉ[2,3] = ν
    
    Cᵉ[3,1] = ν
    Cᵉ[3,2] = ν
    Cᵉ[3,3] = 1 - ν
    
    Cᵉ[4,4] = (1 - 2ν)/2
    
    Cᵉ[5,5] = (1 - 2ν)/2
    
    Cᵉ[6,6] = (1 - 2ν)/2
    
    Cᵉ = E/((1 + ν)*(1 - 2ν)) * Cᵉ
    
    # dielectric tensor
    ϵₜ = zeros(3,3)
    
    ϵₜ[1,1] = ϵ₁
    ϵₜ[2,2] = ϵ₁
    ϵₜ[3,3] = ϵ₃
    
    # magnetic tensor
    μⁱ = zeros(3,3)
    
    μⁱ[1,1] = μᶜ
    μⁱ[2,2] = μᶜ
    μⁱ[3,3] = μᶜ
    
    # piezoelectric tensor
    eₚ = zeros(3,6)
    
    eₚ[1,5] = eₚ14
    eₚ[2,6] = -eₚ14 
    
    return CorticalBone(E, ν, Cᵉ, ϵ₁, ϵ₃, ϵₜ, μᶜ, μⁱ, eₚ14, eₚ)
end


# --- bone marrow
struct BoneMarrow{T, S <: AbstractArray{T, 2}} # T = real (number type)
    E::T # Youngs modulus
    ν::T # Poisson ratio
    Cᵉ::S # Elastic stiffness tensor
    ϵ₁::T # dielectric constant 11
    ϵ₃::T # dielectric constant 33
    ϵₜ::S #  dielectric tensor
    μᶜ::T # inverse of magnetic constant
    μⁱ::S # inverse of magnetic tensor
    κ₁::T # conductivity constant
    κ::S # conductivity tensor
    μᵥ::T # viscoelastic damping parameter (=delta_t*r_1~)
end

function BoneMarrow(E::Float64, ν::Float64, ϵ₁::Float64, ϵ₃::Float64, μᶜ::Float64, κ₁::Float64, μᵥ::Float64) # Constructor zum struct
    # elastic tensor
    Cᵉ = zeros(6,6)
    
    Cᵉ[1,1] = 1 - ν
    Cᵉ[1,2] = ν
    Cᵉ[1,3] = ν
    
    Cᵉ[2,1] = ν
    Cᵉ[2,2] = 1 - ν
    Cᵉ[2,3] = ν
    
    Cᵉ[3,1] = ν
    Cᵉ[3,2] = ν
    Cᵉ[3,3] = 1 - ν
    
    Cᵉ[4,4] = (1 - 2ν)/2
    
    Cᵉ[5,5] = (1 - 2ν)/2
    
    Cᵉ[6,6] = (1 - 2ν)/2
    
    Cᵉ = E/((1 + ν)*(1 - 2ν)) * Cᵉ
    
    # dielectric tensor
    ϵₜ = zeros(3,3)
    
    ϵₜ[1,1] = ϵ₁
    ϵₜ[2,2] = ϵ₁
    ϵₜ[3,3] = ϵ₃
    
    # magnetic tensor
    μⁱ = zeros(3,3)
    
    μⁱ[1,1] = μᶜ
    μⁱ[2,2] = μᶜ
    μⁱ[3,3] = μᶜ
    
    # conductivity tensor
    κ = zeros(3,3)
    
    κ[1,1] = κ₁
    κ[2,2] = κ₁
    κ[3,3] = κ₁
    
    return BoneMarrow(E, ν, Cᵉ, ϵ₁, ϵ₃, ϵₜ, μᶜ, μⁱ, κ₁, κ, μᵥ)
end


# --- state variables
mutable struct MaterialState{T, S <: AbstractArray{T, 1}} # Ergänzen: Felder für die numerischen Störungen 6xsigma,2x3xD,3xJ,3xH
    # Store "converged" values of stress and inelastic strain
    σ::S # stress
    εⁱ::S # inelastic strain
    
    # Store temporary values used during equilibrium iterations (only necessary for σ, εⁱ)
    temp_σ::S
    temp_εⁱ::S
    
    # other state variables
    ε::S # strain
    
    E::S # electric field strength 
    D::S # electric displacement field
    Dp::S # time derivative electric displacement field
    
    B::S # magnetic flux density
    H::S # magnetic field strength
    
    J::S # electric current density
end

function MaterialState()
    return MaterialState(zeros(6),zeros(6),zeros(6),zeros(6),zeros(6),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3))
end

function update_state!(state::MaterialState)
    state.σ = state.temp_σ
    state.εⁱ = state.temp_εⁱ
end

#--- simulation parameters
struct MicroProblem
    a::Float64 # RVE length
    b::Float64 # RVE length
    n::Int64 # RVE mesh resolution
    msp::MacroSimulationParameters
    ctx::CoherentStructures.GridContext # RVE-mesh and utilities <- kann auch nochmal problematisch werden!
    bd::BoundaryData
    S_red::SparseMatrixCSC{Float64,Int64} # stiffness matrix RVE
    nodedoflist::Array{Int64,2}
end

function MicroProblem(micromodel::String, a::Float64, b::Float64, n, msp::MacroSimulationParameters, mp_b::CorticalBone, mp_m::BoneMarrow, etype_micro::DataType) # nur einmal für alle RVEs!
    
    path = string("mesh/",micromodel,".msh")
    grid = read_msh_file(path)
    
    addfaceset!(grid, "left", x -> x[1] ≈ -a-b/2)
    addfaceset!(grid, "right", x -> x[1] ≈ a+b/2)
    addfaceset!(grid, "top", x -> x[2] ≈ -a-b/2)
    addfaceset!(grid, "bottom", x -> x[2] ≈ a+b/2)
    addfaceset!(grid, "front", x -> x[3] ≈ -a-b/2)
    addfaceset!(grid, "back", x -> x[3] ≈ a+b/2)
    addnodeset!(grid, "corners", x -> norm(x[1]) ≈ a+b/2 && norm(x[2]) ≈ a+b/2 && norm(x[3]) ≈ a+b/2)

    ll = Vec{3}((-a-b/2, -a-b/2,-a-b/2))
    ur = Vec{3}((a+b/2, a+b/2, a+b/2))
    loc = CoherentStructures.Regular3DGridLocator{etype_micro}((n+1), (n+1), (n+1), ll, ur)

    dh2 = create_dofhandler_helper(grid, etype_micro)
    
    n_npe = 0 # number of nodes per element
    if(etype_micro == Hexahedron)
        ctx = CoherentStructures.GridContext{3}(grid, Lagrange{3,RefCube,1}(), Lagrange{3,RefCube,1}(), dh2, QuadratureRule{3,RefCube}(2), loc)
        n_npe = 8
    elseif(etype_micro == Tetrahedron)
        ctx = CoherentStructures.GridContext{3}(grid, Lagrange{3,RefTetrahedron,1}(), Lagrange{3,RefTetrahedron,1}(), dh2, QuadratureRule{3,RefTetrahedron}(2), loc)
        n_npe = 4 
    else
        println("Error: micro element type not supported")
    end
        
    grid = nothing
    
    predicate = (p1, p2) -> peuclidean(p1, p2, [2*a+b, 2*a+b, 2*a+b]) < 1e-8 
    bd = BoundaryData(ctx, predicate)
    
    dh = create_dofhandler(ctx.grid, etype_micro)
    cellvalues_u, cellvalues_phi, cellvalues_A = create_values(etype_micro)

    dbc = create_bc(dh)
    S = create_sparsity_pattern(dh)
    
    nodedoflist = nodedofs(dh, ctx.grid)
    
    S = doassemble_S(cellvalues_u, cellvalues_phi, cellvalues_A, S, ctx.grid, dh, mp_b, mp_m, msp)
    n_dofs = ndofs(dh)
    n_dofs_n = Int(ndofs_per_cell(dh)/n_npe) # DoF per node
    r = zeros(n_dofs) # helper vector
    apply_zero!(S, r, dbc) 
    S_red = applyPBC_S(ctx, bd, S, nodedoflist, n_dofs_n)

    return MicroProblem(a, b, n, msp, ctx, bd, S_red, nodedoflist)
end

struct MacroQuantities{T, S <: AbstractArray{T, 1}}
    ε_macro::S # macro strain
    E_macro::S # macro electric field
    B_macro::S # macro magnetic flux density
    t::Int64 # macro time
end