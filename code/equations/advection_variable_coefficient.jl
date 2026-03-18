using Trixi
using Trixi: AbstractCompressibleEulerEquations, norm, stolarsky_mean
import Trixi:
    varnames,
    cons2cons,
    cons2prim,
    boundary_condition_slip_wall,
    have_nonconservative_terms,
    prim2cons,
    cons2entropy, max_abs_speeds, entropy, energy_total

# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

struct AdvectionVariableCoefficientEquation1D{RealT<:Real} <:
       AbstractCompressibleEulerEquations{1,2}
    dummy::RealT
end

function AdvectionVariableCoefficientEquation1D()
    return AdvectionVariableCoefficientEquation1D(0)
end

Trixi.have_nonconservative_terms(::AdvectionVariableCoefficientEquation1D) = Trixi.True()

varnames(::typeof(cons2cons), ::AdvectionVariableCoefficientEquation1D) = ("scalar", "coefficient")
varnames(::typeof(cons2prim), ::AdvectionVariableCoefficientEquation1D) = ("scalar", "coefficient")

@inline function flux_entropy_conservative(u_ll, u_rr, orientation, equations::AdvectionVariableCoefficientEquation1D)
    v1_ll, a_ll = u_ll
    v1_rr, a_rr = u_rr
    jump_u = v1_rr - v1_ll

    flux_left = a_ll * jump_u
    zeroT = zero(eltype(u_ll))
    return SVector(flux_left, zeroT)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Int,
    equation::AdvectionVariableCoefficientEquation1D)
    return abs(equation.advection_velocity[orientation])
end

@inline function max_abs_speeds(equation::AdvectionVariableCoefficientEquation1D)
    return abs.(equation.advection_velocity)
end

@inline function max_abs_speeds(u, equation::AdvectionVariableCoefficientEquation1D)
    return abs.(u[2])
end

# Convert conservative variables to primitive
@inline cons2prim(u, equation::AdvectionVariableCoefficientEquation1D) = u

# Convert conservative variables to entropy variables
@inline cons2entropy(u, equation::AdvectionVariableCoefficientEquation1D) = SVector(2 * u[1] / u[2], zero(eltype(u)))

# Calculate entropy for a conservative state `cons`
@inline function entropy(u, ::AdvectionVariableCoefficientEquation1D)
    return u[1]^2 / u[2]
end

# Calculate total energy for a conservative state `cons`
@inline energy_total(u::Real, ::AdvectionVariableCoefficientEquation1D) = 0.5f0 * u^2
@inline function energy_total(u, equation::AdvectionVariableCoefficientEquation1D)
    energy_total(u[1], equation)
end
end # @muladd
