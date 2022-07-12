# --- macro quantities
struct MacroTangentModuli{T, S <: AbstractArray{T, 2}} # T = real (number type)
    C::S # Mechanic stiffness tensor
    ϵₜ::S #  dielectric tensor
    μⁱ::S # inverse of magnetic tensor 
    eₚ::S # piezoelectric tensor
    κ::S # conductivity tensor
end

# --- macro state variables
mutable struct MacroMaterialState{T, S <: AbstractArray{T, 1}}
    σ::S # stress
    ε::S # strain
    E::S # electric field strength 
    D::S # electric displacement field
    Dp::S # time derivative electric displacement field
    B::S # magnetic flux density
    H::S # magnetic field strength
    J::S # electric current density
end

function MacroMaterialState()
    return MacroMaterialState(zeros(6),zeros(6),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3))
end

#--- macro simulation parameters
struct MacroSimulationParameters{T, S <: AbstractArray{T, 1}} 
    # time integration parameters
    ρ_∞::T 
    α_f::T
    αₘ::T
    γₐ::T

    Δₜ::T # time increment
    NEWTON_TOL::T # tolerance for convergence check of NR-method
    γ::T # punish parameter for divergence of magnetic vector potential
    ctan::S # weight vector for time integration
end

function MacroSimulationParameters(ρ_∞::Float64, Δₜ::Float64, NEWTON_TOL::Float64, γ::Float64)
    # timesteps / time integration (Gl.32) - JWH-alpha method Kadapa 2017
    # Verfahrensparameter (Gl. 32)
    # ρ_∞ = 0.5 # 0 <= rho_inf <= 1, ggbf eher klein halten, Rest abh. davon
    α_f = 1/(1+ρ_∞)
    αₘ = (3-ρ_∞)/(2*(1+ρ_∞))
    γₐ = 0.5 + αₘ - α_f
    ctan = [α_f, αₘ/(γₐ*Δₜ), αₘ^2/(α_f*γₐ^2*Δₜ^2)]
    return MacroSimulationParameters(ρ_∞, α_f, αₘ, γₐ, Δₜ, NEWTON_TOL, γ, ctan)
end

