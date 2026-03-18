using Trixi
using Trixi: AbstractCompressibleEulerEquations, norm, stolarsky_mean
import Trixi:
    varnames,
    cons2cons,
    cons2prim,
    boundary_condition_slip_wall,
    have_nonconservative_terms,
    prim2cons,
    cons2entropy, max_abs_speeds, entropy

# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
    PolynomialEquation1D

The linear scalar advection with varible coefficient equation
```math
\partial_t u + a(x) \partial_1 u  = 0
```
in one space dimension with velocity `a(x)`.
"""
struct PolynomialEquation1D{RealT <: Real} <: AbstractCompressibleEulerEquations{1,1}
    m::Int
    n::Int
    alpha::RealT
end

function PolynomialEquation1D(m, n , alpha)
    PolynomialEquation1D(m, n, alpha)
end

Trixi.have_nonconservative_terms(::PolynomialEquation1D) = Trixi.True()

@inline function flux_zero(u_ll, u_rr, orientation, equations)
    zeroT = zero(eltype(u_ll))
    return SVector(zeroT)
end

@inline function flux_entropy_conservative(
    u_ll,
    u_rr,
    orientation::Integer,
    equations::PolynomialEquation1D,
)
    u_ll = u_ll[1]
    u_rr = u_rr[1]

    jump = u_rr - u_ll
    alpha = equations.alpha
    k = equations.k
    h_num = 0
    for i in 0:equations.k-1
        h_num += u_rr^(k-1-i) * u_ll^i 
    end
    h_num = h_num/k

    flux_left = alpha * k * h_num * jump * 0.5f0 + (1-alpha) * k * u_ll^(k-1) * jump * 0.5f0

    return SVector(flux_left)
end

@inline function flux_entropy_general(
    u_ll,
    u_rr,
    orientation::Integer,
    equations::PolynomialEquation1D,
)
    u_ll = u_ll[1]
    u_rr = u_rr[1]
    m = equations.m
    n = equations.n
    jump = u_rr^n - u_ll^n
    alpha = equations.alpha
    h_num_jump1 = 0
    for k in 0:(m+n)
        h_num_jump1 += (-1)^k * u_rr^(m+n-k) * u_ll^k 
    end
    h_num_jump2 = 0
    for k in 0:(n-1)
        h_num_jump2 += (-1)^k * u_rr^(n-1-k) * u_ll^k 
    end
    h_num_jump2 = -h_num_jump2 *(u_rr^(m+1) + u_ll^(m+1)) * 0.5f0

    h_num_jump2 = 0
    for k in 0:m
	h_num_jump2 += (-1)^k * u_rr^(m-k) * u_ll^(n+k) 
    end
    h_num_jump3 = 0
    for k in 0:(n-1)
	h_num_jump3 += (-1)^k * u_rr^(n-1-k) * u_ll^(m+1+k) 
    end
    h_num_jump = h_num_jump1 + h_num_jump2 - h_num_jump3
    h_num_jump = h_num_jump * n/(m+1)

    flux_left = alpha * h_num_jump * 0.5f0 + (1-alpha) * u_ll^m * jump * 0.5f0

    return SVector(flux_left)
end

@inline function flux_nonconservative_even(u_ll, u_rr, orientation::Integer, equations::PolynomialEquation1D)

	u_ll = u_ll[1]
	u_rr = u_rr[1]

	m = equations.m
	n = equations.n
   	
    alpha = equations.alpha
    jump = u_rr^m - u_ll^m
    h_num = (u_ll^n + u_rr^n)*0.5f0	
    flux_left = alpha * h_num * jump + (1-alpha) * u_ll^n * jump 
       return SVector(-flux_left)
end


@inline function flux_conservative_even(u_ll, u_rr, orientation::Integer, equations::PolynomialEquation1D)

	u_ll = u_ll[1]
	u_rr = u_rr[1]

	m = equations.m
	n = equations.n
        flux_1 = 0

	for j in 0:(m-1)
	flux_1 =flux_1 + u_rr^(m-1-j)*u_ll^j
	end	
	u_avg = 0.5f0 * (u_ll + u_rr)
	unp1_avg = 0.5f0 * (u_ll^(n+1) + u_rr^(n+1))
	un_avg = 0.5f0 * (u_ll^(n) + u_rr^(n))
	flux_1 = flux_1 * (equations.alpha * un_avg * u_avg + (1-equations.alpha) * unp1_avg)

        flux_2 = 0

	for j in 0:(m+n)
	flux_2 = flux_2 + u_rr^(m+n-j)*u_ll^j
	end	
	flux_2 = flux_2 * (m+1)/(m+n+1)
        flux = -flux_1 + flux_2	
	return SVector(flux)

end


varnames(::typeof(cons2cons), ::PolynomialEquation1D) = ("u")
varnames(::typeof(cons2prim), ::PolynomialEquation1D) = ("u")

# Convert conservative variables to primitive
@inline cons2prim(u, equation::PolynomialEquation1D) = u

# Convert conservative variables to entropy variables
@inline cons2entropy(u, equation::PolynomialEquation1D) = u

# Calculate entropy for a conservative state `cons`
@inline entropy(u::Real, ::PolynomialEquation1D) = 0.5f0 * u^2
@inline entropy(u, equation::PolynomialEquation1D) = entropy(u[1], equation)

# Calculate total energy for a conservative state `cons`
@inline energy_total(u::Real, ::PolynomialEquation1D) = 0.5f0 * u^2
@inline function energy_total(u, equation::PolynomialEquation1D)
    energy_total(u[1], equation)
end

@inline function velocity(u, equation::PolynomialEquation1D)
	return u[1]
end

@inline function Trixi.max_abs_speeds(u, equation::PolynomialEquation1D)
    return (abs(u[1]),)
end
end # @muladd
