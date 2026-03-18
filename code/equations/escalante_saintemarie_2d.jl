@muladd begin
#! format: noindent

struct EscalanteSainteMarieEquations2D{RealT<:Real} <:
       Trixi.AbstractShallowWaterEquations{2,6}
    gravity::RealT
    H0::RealT
    celerity::RealT
    threshold_limiter::RealT
    threshold_wet::RealT
    threshold_partially_wet::RealT
    threshold_desingularization::RealT
end

function EscalanteSainteMarieEquations2D(; gravity, H0=zero(gravity), b0=one(gravity),
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

    EscalanteSainteMarieEquations2D(gravity, H0, celerity, threshold_limiter,
        threshold_wet, threshold_partially_wet,
        threshold_desingularization)
end

Trixi.have_nonconservative_terms(::EscalanteSainteMarieEquations2D) = Trixi.True()

function Trixi.varnames(::typeof(cons2cons), ::EscalanteSainteMarieEquations2D)
    ("h", "h_u", "h_v", "h_w", "h_p", "b")
end
# Note, we use the total water height, H = h + b, as the first primitive variable for easier
# visualization and setting initial conditions
function Trixi.varnames(::typeof(cons2prim), ::EscalanteSainteMarieEquations2D)
    ("H", "u","v", "w", "p", "b")
end

@inline function Trixi.boundary_condition_slip_wall(u_inner,
                                                    normal_direction::AbstractVector,
                                                    x, t,
                                                    surface_flux_functions,
                                                    equations::EscalanteSainteMarieEquations2D)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions

    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4], u_inner[5], u_inner[6])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)
    return flux, noncons_flux
end

@inline function source_terms_escalante(u, x, t, equations::EscalanteSainteMarieEquations2D)
    du1 = zero(eltype(u))
    du2 = zero(eltype(u))
    du3 = zero(eltype(u))
    du4 = 2 * u[5] / u[1]
    du5 = -2 * equations.celerity^2 * u[4] / u[1]
    du6 = zero(eltype(u))

    return SVector(du1, du2, du3, du4, du5, du6)
end

@inline function flux_conservative(u_ll, u_rr,
    normal_direction::AbstractVector,
    equations::EscalanteSainteMarieEquations2D)
    # Pull the necessary left and right state information
    h_ll, h_v1_ll, h_v2_ll, h_v3_ll, h_p_ll, b_ll = u_ll
    h_rr, h_v1_rr, h_v2_rr, h_v3_rr, h_p_rr, b_rr = u_rr

    v1_ll = h_v1_ll / h_ll
    v2_ll = h_v2_ll / h_ll
    v3_ll = h_v3_ll / h_ll
    p_ll = h_p_ll / h_ll

    v1_rr = h_v1_rr / h_rr
    v2_rr = h_v2_rr / h_rr
    v3_rr = h_v3_rr / h_rr
    p_rr = h_p_rr / h_rr

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    v3_avg = 0.5f0 * (v3_ll + v3_rr)

    p_avg = 0.5f0 * (p_ll + p_rr)

    h_avg = 0.5f0 * (h_ll + h_rr)
    h_v1_avg = 0.5f0 * (h_v1_ll + h_v1_rr)
    h_v2_avg = 0.5f0 * (h_v2_ll + h_v2_rr)
    h_v_dot_avg = h_v1_avg * normal_direction[1] + h_v2_avg * normal_direction[2]
    h_p_avg = 0.5f0 * (h_p_ll + h_p_rr)
    h2_avg = 0.5f0 * (h_ll^2 + h_rr^2)
 v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f1 = alpha1 * h_avg * v_dot_n_avg + (1 - alpha1) * h_v_dot_avg
    pressure_terms = equations.gravity * (1 - alpha1) * h_avg^2 + equations.gravity * (alpha1 - 0.5f0) * h2_avg + alpha3 * h_avg * p_avg + (1 - alpha3) * h_p_avg
    f2 = f1 * v1_avg + pressure_terms * normal_direction[1]
    f3 = f1 * v2_avg + pressure_terms * normal_direction[2]
    f4 = f1 * v3_avg
    f5 = f1 * p_avg

    f = SVector(f1,
        f2,
        f3, f4, f5, zero(eltype(u_ll)))

    return f
end

@inline function flux_nonconservative(u_ll, u_rr,
    normal_direction::AbstractVector,
    equations::EscalanteSainteMarieEquations2D)
    # Pull the necessary left and right state information
    h_ll, h_v1_ll, h_v2_ll, h_v3_ll, h_p_ll, b_ll = u_ll
    h_rr, h_v1_rr, h_v2_rr, h_v3_rr, h_p_rr, b_rr = u_rr

    v1_ll = h_v1_ll / h_ll
    v2_ll = h_v2_ll / h_ll
    p_ll = h_p_ll / h_ll
    b_jump = b_rr - b_ll
    v1_rr = h_v1_rr / h_rr
    v2_rr = h_v2_rr / h_rr
    p_rr = h_p_rr / h_rr

    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    p_avg = 0.5f0 * (p_ll + p_rr)

    h_avg = 0.5f0 * (h_ll + h_rr)

    v1_jump = v1_rr - v1_ll
    v2_jump = v2_rr - v2_ll

    f1 = zero(eltype(u_ll))
    f4 = zero(eltype(u_ll))
    f2_fluxdiff = alpha1 * equations.gravity * h_avg * b_jump + alpha2 * 2 * p_avg * b_jump
    f2_pointwise = (1 - alpha1) * equations.gravity * h_ll * b_jump + (1 - alpha2) * 2 * p_ll * b_jump
    f2 = f2_fluxdiff + f2_pointwise
    f3 = f2

    f5_fluxdiff = alpha3 * equations.celerity^2 * h_avg * v1_jump - 2 * alpha2 * equations.celerity^2 * v1_avg * b_jump
    f5_pointwise = (1 - alpha3) * equations.celerity^2 * h_ll * v1_jump - 2 * (1 - alpha2) * equations.celerity^2 * v1_ll * b_jump
    f5_v1 = f5_fluxdiff + f5_pointwise

    f5_fluxdiff = alpha3 * equations.celerity^2 * h_avg * v2_jump - 2 * alpha2 * equations.celerity^2 * v2_avg * b_jump
    f5_pointwise = (1 - alpha3) * equations.celerity^2 * h_ll * v2_jump - 2 * (1 - alpha2) * equations.celerity^2 * v2_ll * b_jump
    f5_v2 = f5_fluxdiff + f5_pointwise
    f5 = f5_v1 * normal_direction[1] + f5_v2 * normal_direction[2]
    f = SVector(f1,
        f2 * normal_direction[1],
        f3 * normal_direction[2],
        f4,
        f5,
        zero(eltype(u_ll)))

    return f
end


@inline function Trixi.max_abs_speeds(u, equations::EscalanteSainteMarieEquations2D)
    h = u[1]
    v1, v2 = velocity(u, equations)
    p = u[5] / h
    c = sqrt(equations.gravity * h + p + equations.celerity^2)
    return abs(v1) + c, abs(v2) + c
end

# Helper function to extract the velocity vector from the conservative variables
@inline function Trixi.velocity(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2 = u

    v1 = h_v1 / h
    v2 = h_v2 / h

    return v1, v2
end

# Convert conservative variables to primitive
@inline function Trixi.cons2prim(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, h_p, b = u

    H = h + b
    v1 = h_v1 / h
    v2 = h_v2 / h
    v3 = h_v3 / h
    p = h_p / h
    return SVector(H, v1, v2, v3, p, b)
end

# Convert conservative variables to entropy
# Note, only the first two are the entropy variables, the third entry still
# just carries the bottom topography values for convenience
@inline function Trixi.cons2entropy(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, h_p, b = u

    v1 = h_v1 / h
    v2 = h_v2 / h
    v3 = h_v3 / h
    p = h_p / h

    w1 = equations.gravity * (h + b) - 0.5f0 * (v1^2 + v2^2 + v3^2 + p^2 / equations.celerity^2)
    w2 = v1
    w3 = v2
    w4 = v3
    w5 = p / equations.celerity^2
    return SVector(w1, w2, w3, w4, w5, zero(eltype(u)))
end

# Convert primitive to conservative variables
@inline function Trixi.prim2cons(prim, equations::EscalanteSainteMarieEquations2D)
    H, v1, v2, v3, p, b = prim

    h = H - b
    h_v1 = h * v1
    h_v2 = h * v2
    h_v3 = h * v3
    h_p = h * p

    return SVector(h, h_v1, h_v2, h_v3, h_p, b)
end

# Entropy function for the shallow water equations is the total energy
@inline function Trixi.entropy(cons, equations::EscalanteSainteMarieEquations2D)
    energy_total(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline function Trixi.energy_total(cons, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, h_p, b = cons

    e = (h_v1^2 + h_v2^2 + h_v3^2) / (2 * h) + h_p^2 / (2 * equations.celerity^2 * h) + 0.5f0 * equations.gravity * h^2 + equations.gravity * h * b
    return e
end

@inline function Trixi.lake_at_rest_error(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, hp, b = u
    H0_wet_dry = max(equations.H0, b + equations.threshold_limiter)

    return abs(H0_wet_dry - (h + b))
end

@inline function velocity_1(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, hp, b = u

    return h_v1/h
end

@inline function velocity_2(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, hp, b = u

    return h_v2/h
end

@inline function velocity_3(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, hp, b = u

    return h_v3/h
end

@inline function pressure_(u, equations::EscalanteSainteMarieEquations2D)
    h, h_v1, h_v2, h_v3, h_p, b = u

    return h_p/h
end

end # @muladd
