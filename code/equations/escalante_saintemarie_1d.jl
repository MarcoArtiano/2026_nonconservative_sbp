@muladd begin
#! format: noindent

struct EscalanteSainteMarieEquations1D{RealT<:Real} <:
       Trixi.AbstractShallowWaterEquations{1,5}
    gravity::RealT
    H0::RealT
    celerity::RealT
    threshold_limiter::RealT
    threshold_wet::RealT
    threshold_partially_wet::RealT
    threshold_desingularization::RealT
end

function EscalanteSainteMarieEquations1D(; gravity, H0=zero(gravity), b0=one(gravity),
    threshold_limiter=nothing,
    threshold_wet=nothing,
    threshold_partially_wet=nothing,
    threshold_desingularization=nothing)
    T = promote_type(typeof(gravity), typeof(H0), typeof(b0))
    if threshold_limiter === nothing
        threshold_limiter = 500 * eps(T)
    end
    if threshold_wet === nothing
        threshold_wet = 5 * eps(T)
    end
    if threshold_partially_wet === nothing
        threshold_partially_wet = TrixiShallowWater.default_threshold_partially_wet(T)
    end
    if threshold_desingularization === nothing
        threshold_desingularization = TrixiShallowWater.default_threshold_desingularization(T)
    end
    celerity = 2 * sqrt(gravity * b0)

    EscalanteSainteMarieEquations1D(gravity, H0, celerity, threshold_limiter,
        threshold_wet, threshold_partially_wet,
        threshold_desingularization)
end

Trixi.have_nonconservative_terms(::EscalanteSainteMarieEquations1D) = Trixi.True()

function Trixi.varnames(::typeof(cons2cons), ::EscalanteSainteMarieEquations1D)
    ("h", "h_v", "h_w", "h_p", "b")
end
# Note, we use the total water height, H = h + b, as the first primitive variable for easier
# visualization and setting initial conditions
function Trixi.varnames(::typeof(cons2prim), ::EscalanteSainteMarieEquations1D)
    ("H", "v", "w", "p", "b")
end

"""
    boundary_condition_slip_wall(u_inner, orientation_or_normal, x, t, surface_flux_function,
                                  equations::EscalanteSainteMarieEquations1D)

Create a boundary state by reflecting the normal velocity component and keep
the tangential velocity component unchanged. The boundary water height is taken from
the internal value.

For details see Section 9.2.5 of the book:
- Eleuterio F. Toro (2001)
  Shock-Capturing Methods for Free-Surface Shallow Flows
  1st edition
  ISBN 0471987662
"""
@inline function Trixi.boundary_condition_slip_wall(u_inner, orientation_or_normal,
    direction,
    x, t,
    surface_flux_functions,
    equations::EscalanteSainteMarieEquations1D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # This can not be dispatched, when Flux Hydrostactic reconstruction is used
    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
        -u_inner[2],
        u_inner[3], u_inner[4], u_inner[5])

    # calculate the boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal,
            equations)
        noncons_flux = nonconservative_flux_function(u_inner, u_boundary,
            orientation_or_normal,
            equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal,
            equations)
        noncons_flux = nonconservative_flux_function(u_boundary, u_inner,
            orientation_or_normal,
            equations)
    end

    return flux, noncons_flux
end

@inline function source_terms_escalante(u, x, t, equations::EscalanteSainteMarieEquations1D)
    du1 = zero(eltype(u))
    du2 = zero(eltype(u))
    du3 = 2 * u[4] / u[1]
    du4 = -2 * equations.celerity^2 * u[3] / u[1]
    du5 = zero(eltype(u))

    return SVector(du1, du2, du3, du4, du5)
end

@inline function flux_conservative(u_ll, u_rr,
    orientation::Integer,
    equations::EscalanteSainteMarieEquations1D)
    # Pull the necessary left and right state information
    h_ll, h_v1_ll, h_v2_ll, h_p_ll, b_ll = u_ll
    h_rr, h_v1_rr, h_v2_rr, h_p_rr, b_rr = u_rr

    v1_ll = h_v1_ll / h_ll
    v2_ll = h_v2_ll / h_ll
    p_ll = h_p_ll / h_ll

    v1_rr = h_v1_rr / h_rr
    v2_rr = h_v2_rr / h_rr
    p_rr = h_p_rr / h_rr

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    h_avg = 0.5f0 * (h_ll + h_rr)
    h_v1_avg = 0.5f0 * (h_v1_ll + h_v1_rr)
    h_p_avg = 0.5f0 * (h_p_ll + h_p_rr)
    h2_avg = 0.5f0 * (h_ll^2 + h_rr^2)

    f1 = alpha1 * h_avg * v1_avg + (1 - alpha1) * h_v1_avg
    pressure_terms = equations.gravity * (1 - alpha1) * h_avg^2 + equations.gravity * (alpha1 - 0.5f0) * h2_avg + alpha3 * h_avg * p_avg + (1 - alpha3) * h_p_avg
    f2 = f1 * v1_avg + pressure_terms
    f3 = f1 * v2_avg
    f4 = f1 * p_avg

    f = SVector(f1,
        f2,
        f3, f4, zero(eltype(u_ll)))

    return f
end

@inline function flux_nonconservative(u_ll, u_rr,
    orientation::Integer,
    equations::EscalanteSainteMarieEquations1D)
    # Pull the necessary left and right state information
    h_ll, h_v1_ll, h_v2_ll, h_p_ll, b_ll = u_ll
    h_rr, h_v1_rr, h_v2_rr, h_p_rr, b_rr = u_rr

    v1_ll = h_v1_ll / h_ll
    v2_ll = h_v2_ll / h_ll
    p_ll = h_p_ll / h_ll
    b_jump = b_rr - b_ll
    v1_rr = h_v1_rr / h_rr
    p_rr = h_p_rr / h_rr

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    h_avg = 0.5f0 * (h_ll + h_rr)

    v1_jump = v1_rr - v1_ll

    f1 = zero(eltype(u_ll))
    f2_fluxdiff = alpha1 * equations.gravity * h_avg * b_jump + alpha2 * 2 * p_avg * b_jump
    f2_pointwise = (1 - alpha1) * equations.gravity * h_ll * b_jump + (1 - alpha2) * 2 * p_ll * b_jump
    f2 = f2_fluxdiff + f2_pointwise
    f3 = zero(eltype(u_ll))

    f4_fluxdiff = alpha3 * equations.celerity^2 * h_avg * v1_jump - 2 * alpha2 * equations.celerity^2 * v1_avg * b_jump
    f4_pointwise = (1 - alpha3) * equations.celerity^2 * h_ll * v1_jump - 2 * (1 - alpha2) * equations.celerity^2 * v1_ll * b_jump
    f4 = f4_fluxdiff + f4_pointwise

    f = SVector(f1,
        f2,
        f3,
        f4,
        zero(eltype(u_ll)))

    return f
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function Trixi.max_abs_speed(u_ll, u_rr, orientation::Integer,
    equations::EscalanteSainteMarieEquations1D)
    # Get the velocity quantities
    v_ll = velocity(u_ll, equations)
    v_rr = velocity(u_rr, equations)

    # Calculate the wave celerity on the left and right
    h_ll = u_ll[1]
    h_rr = u_rr[1]
    p_ll = u_ll[4] / h_ll
    p_rr = u_rr[4] / h_rr
    c_ll = sqrt(equations.gravity * h_ll + p_ll + equations.celerity^2)
    c_rr = sqrt(equations.gravity * h_rr + p_rr + equations.celerity^2)

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end


@inline function Trixi.max_abs_speeds(u, equations::EscalanteSainteMarieEquations1D)
    h = u[1]
    v = velocity(u, equations)
    p = u[4] / h
    c = sqrt(equations.gravity * h + p + equations.celerity^2)
    return (abs(v) + c,)
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::EscalanteSainteMarieEquations1D)
    h, h_v, _ = u

    v = h_v / h

    return v
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::EscalanteSainteMarieEquations1D)
    h, h_v1, h_v2, h_p, b = u

    H = h + b
    v1 = h_v1 / h
    v2 = h_v2 / h
    p = h_p / h
    return SVector(H, v1, v2, p, b)
end

# Convert conservative variables to entropy
# Note, only the first two are the entropy variables, the third entry still
# just carries the bottom topography values for convenience
@inline function Trixi.cons2entropy(u, equations::EscalanteSainteMarieEquations1D)
    h, h_v1, h_v2, h_p, b = u

    v1 = h_v1 / h
    v2 = h_v2 / h
    p = h_p / h

    w1 = equations.gravity * (h + b) - 0.5f0 * (v1^2 + v2^2 + p^2 / equations.celerity^2)
    w2 = v1
    w3 = v2
    w4 = p / equations.celerity^2
    return SVector(w1, w2, w3, w4, zero(eltype(u)))
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::EscalanteSainteMarieEquations1D)
    H, v1, v2, p, b = prim

    h = H - b
    h_v1 = h * v1
    h_v2 = h * v2
    h_p = h * p

    return SVector(h, h_v1, h_v2, h_p, b)
end

# Entropy function for the shallow water equations is the total energy
@inline function Trixi.entropy(cons, equations::EscalanteSainteMarieEquations1D)
    energy_total(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline function Trixi.energy_total(cons, equations::EscalanteSainteMarieEquations1D)
    h, h_v1, h_v2, h_p, b = cons

    e = (h_v1^2 + h_v2^2) / (2 * h) + h_p^2 / (2 * equations.celerity^2 * h) + 0.5f0 * equations.gravity * h^2 + equations.gravity * h * b
    return e
end

@inline function Trixi.lake_at_rest_error(u, equations::EscalanteSainteMarieEquations1D)
    h, _, _, _, b = u
    H0_wet_dry = max(equations.H0, b + equations.threshold_limiter)

    return abs(H0_wet_dry - (h + b))
end
end # @muladd
