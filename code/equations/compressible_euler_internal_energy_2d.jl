using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations, norm
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, FluxLMARS, boundary_condition_slip_wall, have_nonconservative_terms, prim2cons, energy_total

@muladd begin
#! format: noindent
struct CompressibleEulerInternalEnergyEquationsWithGravity2D{RealT <: Real} <:
       AbstractCompressibleEulerEquations{2, 5}
    p_0::RealT # reference pressure in Pa
    c_p::RealT # specific heat at constant pressure in J/(kg K)
    c_v::RealT # specific heat at constant volume in J/(kg K)
    g::RealT # gravitational acceleration in m/s²
    R::RealT # gas constant in J/(kg K)
    gamma::RealT # ratio of specific heats
    inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications
    K::RealT # = p_0 * (R / p_0)^gamma; scaling factor between pressure and weighted potential temperature
    stolarsky_factor::RealT # = (gamma - 1) / gamma; used in the stolarsky mean
    function CompressibleEulerInternalEnergyEquationsWithGravity2D(; c_p, c_v,
                                                                         gravity)
        c_p, c_v, g = promote(c_p, c_v, gravity)
        p_0 = 100_000
        R = c_p - c_v
        gamma = c_p / c_v
        inv_gamma_minus_one = inv(gamma - 1)
        K = p_0 * (R / p_0)^gamma
        stolarsky_factor = (gamma - 1) / gamma
        return new{typeof(c_p)}(p_0, c_p, c_v, g, R,
                                gamma,
                                inv_gamma_minus_one,
                                K, stolarsky_factor)
    end
end

function varnames(::typeof(cons2cons),
                  ::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    ("rho", "rho_v1", "rho_v2", "rho_e", "phi")
end

varnames(::typeof(cons2prim),
::CompressibleEulerInternalEnergyEquationsWithGravity2D) = ("rho",
                                                                  "v1",
                                                                  "v2",
                                                                  "p", "phi")

have_nonconservative_terms(::CompressibleEulerInternalEnergyEquationsWithGravity2D) = Trixi.True()

@inline function boundary_condition_slip_wall(u_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
                                              surface_flux_functions,
                                              equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    surface_flux_function, nonconservative_flux_function = surface_flux_functions
    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] 

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4], u_inner[5])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
                                                 equations)
    return flux, noncons_flux
end



@inline function flux_surface_combined(u_ll, u_rr,
                                                       normal_direction::AbstractVector,
                                                       equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    rho_mean = ln_mean(rho_ll, rho_rr)
    phi_jump = phi_rr - phi_ll
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the necessary mean values
end


@inline function flux_surface_combined(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    v_dot_n_avg = 0.5f0 * (v_dot_n_rr + v_dot_n_ll)

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    a = 340.0
    norm_ = norm(normal_direction)
    T_ll = p_ll/(equations.R * rho_ll)
    T_rr = p_rr/(equations.R * rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    p_interface = 0.5f0 * (p_ll + p_rr) -0.5f0 * a * rho_avg * (v_dot_n_rr - v_dot_n_ll)/norm_
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) -1/(2 * a * rho_avg) * (p_rr - p_ll) * norm_
    p_avg = 0.5f0 * (p_ll + p_rr)
    if (v_interface >=  0)
        rho_up = rho_mean - 0.5f0 * (rho_rr - rho_ll) 
        f1, f2, f3, f4 = u_ll * v_interface
    else
        rho_up = rho_mean + 0.5f0 * (rho_rr - rho_ll) 
  f1, f2, f3, f4 = u_rr * v_interface
    end
    T_avg = Trixi.inv_ln_mean(1 / T_ll, 1 / T_rr) * equations.c_v
    # Calculate fluxes depending on normal_direction
    f1 = rho_up * v_interface
   # f2 = f1 * v1_avg + p_interface * normal_direction[1] 
   # f3 = f1 * v2_avg + p_interface * normal_direction[2]
    f2 = f2 + p_interface * normal_direction[1] 
    f3 = f3 + p_interface * normal_direction[2]
    f4 = f1 * T_avg + p_avg * v_interface
    flux = SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
    phi_jump = phi_rr - phi_ll
    gravity = 0.5f0 * rho_mean * phi_jump
    flux_internal_energy = - 0.5f0 * v_dot_n_avg * (p_rr - p_ll)
    noncons_left = SVector(zero(eltype(u_ll)), normal_direction[1] * gravity, normal_direction[2] * gravity, flux_internal_energy, zero(eltype(u_ll)))
    noncons_right = -noncons_left
    flux_left = flux + noncons_left
    flux_right = flux + noncons_right
    return flux_left, flux_right
end


 @inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_surface_combined),
    equations::CompressibleEulerInternalEnergyEquationsWithGravity2D) = Trixi.True()


@inline function boundary_condition_slip_wall(u_inner,
                                              normal_direction::AbstractVector,
                                              x, t,
					      surface_flux_function::typeof(flux_surface_combined),
                                              equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)
    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
                         u_inner[2] - 2 * u_normal * normal[1],
                         u_inner[3] - 2 * u_normal * normal[2],
                         u_inner[4], u_inner[5])

    # calculate the boundary flux
    flux, _ = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
    return flux
end

@inline function flux_volume_ec_combined(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    v_dot_n_avg = 0.5f0 * (v_dot_n_rr + v_dot_n_ll)
    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)

    p_interface = 0.5f0 * (p_ll + p_rr)
    T_ll = p_ll/(equations.R * rho_ll)
    T_rr = p_rr/(equations.R * rho_rr)
    e_avg = Trixi.inv_ln_mean(1 / T_ll, 1 / T_rr) * equations.c_v
    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * v_dot_n_avg
    f2 = f1 * v1_avg + p_interface * normal_direction[1] 
    f3 = f1 * v2_avg + p_interface * normal_direction[2]
    f4 = f1 * e_avg + p_interface * v_dot_n_avg
    flux = SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
    phi_jump = phi_rr - phi_ll
    gravity = 0.5f0 * rho_mean * phi_jump
    flux_internal_energy = - 0.5f0 * v_dot_n_avg * (p_rr - p_ll)
    noncons_left = SVector(zero(eltype(u_ll)), normal_direction[1] * gravity, normal_direction[2] * gravity, flux_internal_energy, zero(eltype(u_ll)))
    noncons_right = -noncons_left
    flux_left = flux + noncons_left
    flux_right = flux + noncons_right
    return flux_left, flux_right
    end

 @inline Trixi.combine_conservative_and_nonconservative_fluxes(::typeof(flux_volume_ec_combined),
    equations::CompressibleEulerInternalEnergyEquationsWithGravity2D) = Trixi.True()

@inline function flux_nonconservative_ec(u_ll, u_rr,
                                                       normal_direction::AbstractVector,
                                                       equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    rho_avg = ln_mean(rho_ll, rho_rr)
    phi_jump = phi_rr - phi_ll
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    p_avg = 0.5f0 * (p_rr + p_ll)
    flux_internal_energy = p_avg * (v_dot_n_rr - v_dot_n_ll)
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   flux_internal_energy, zero(eltype(u_ll)))
end


@inline function flux_ec(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)

    p_interface = 0.5f0 * (p_ll + p_rr)
    T_ll = p_ll/(equations.R * rho_ll)
    T_rr = p_rr/(equations.R * rho_rr)
    T_avg = Trixi.inv_ln_mean(1 / T_ll, 1 / T_rr) * equations.c_v
    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_interface * normal_direction[1] 
    f3 = f1 * v2_avg + p_interface * normal_direction[2]
    f4 = f1 * T_avg 
    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

@inline function flux_nonconservative_ec1(u_ll, u_rr,
                                                       normal_direction::AbstractVector,
                                                       equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    rho_mean = ln_mean(rho_ll, rho_rr) 
    phi_jump = phi_rr - phi_ll
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    a = 340.0
    norm_ = norm(normal_direction)

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) -1/(2 * a * rho_avg) * (p_rr - p_ll) * norm_
    p_avg = 0.5f0 * (p_rr + p_ll)
    flux_internal_energy = -v_interface * (p_rr - p_ll)
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_mean * phi_jump,
                   normal_direction[2] * rho_mean * phi_jump,
                   flux_internal_energy, zero(eltype(u_ll)))
end


@inline function flux_ec1(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the necessary mean values
    rho_mean = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)

    p_interface = 0.5f0 * (p_ll + p_rr)
    T_ll = p_ll/(equations.R * rho_ll)
    T_rr = p_rr/(equations.R * rho_rr)
    T_avg = 0.5f0 * (T_ll + T_rr) * equations.c_v
    p_avg = 0.5f0 * (p_rr + p_ll)
    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_interface * normal_direction[1] 
    f3 = f1 * v2_avg + p_interface * normal_direction[2]
    f4 = f1 * T_avg + p_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) 
    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end

@inline function flux_nonconservative_entropy_stable(u_ll, u_rr,
                                                       normal_direction::AbstractVector,
                                                       equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    rho_mean = ln_mean(rho_ll, rho_rr)
    phi_jump = phi_rr - phi_ll
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the necessary mean values
    a = 340.0
    norm_ = norm(normal_direction)

    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) -1/(2 * a * rho_avg) * (p_rr - p_ll) * norm_

    flux_internal_energy = -v_interface * (p_rr - p_ll)
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_mean * phi_jump,
                   normal_direction[2] * rho_mean * phi_jump,
                   flux_internal_energy, zero(eltype(u_ll)))
end


@inline function flux_entropy_stable(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the necessary mean values
    rho_mean = ln_mean(rho_ll, rho_rr)
    a = 340.0
    norm_ = norm(normal_direction)
    T_ll = p_ll/(equations.R * rho_ll)
    T_rr = p_rr/(equations.R * rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    p_interface = 0.5f0 * (p_ll + p_rr) -0.5f0 * a * rho_avg * (v_dot_n_rr - v_dot_n_ll)/norm_
    v_interface = 0.5f0 * (v_dot_n_ll + v_dot_n_rr) -1/(2 * a * rho_avg) * (p_rr - p_ll) * norm_
    p_avg = 0.5f0 * (p_ll + p_rr)
    if (v_interface >=  0)
        rho_up = rho_mean - 0.5f0 * (rho_rr - rho_ll) 
        f1, f2, f3, f4 = u_ll * v_interface
    else
        rho_up = rho_mean + 0.5f0 * (rho_rr - rho_ll) 
  f1, f2, f3, f4 = u_rr * v_interface
    end
    T_avg = Trixi.inv_ln_mean(1 / T_ll, 1 / T_rr) * equations.c_v
    # Calculate fluxes depending on normal_direction
    f1 = rho_up * v_interface
   # f2 = f1 * v1_avg + p_interface * normal_direction[1] 
   # f3 = f1 * v2_avg + p_interface * normal_direction[2]
    f2 = f2 + p_interface * normal_direction[1] 
    f3 = f3 + p_interface * normal_direction[2]
    f4 = f1 * T_avg + p_avg * v_interface
    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end




@inline function flux_nonconservative_avg(u_ll, u_rr,
                                                       normal_direction::AbstractVector,
                                                       equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
    rho_avg = 0.5f0 * (rho_ll + rho_rr)
    phi_jump = phi_rr - phi_ll
   
    flux_internal_energy = -v_ll * (p_rr - p_ll)
    return SVector(zero(eltype(u_ll)),
                   normal_direction[1] * rho_avg * phi_jump,
                   normal_direction[2] * rho_avg * phi_jump,
                   flux_internal_energy, zero(eltype(u_ll)))
end


@inline function flux_avg(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
    v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    # Compute the necessary mean values
    rho_mean = 0.5f0 * (rho_ll + rho_rr)
    v1_avg = 0.5f0 * (v1_ll + v1_rr)
    v2_avg = 0.5f0 * (v2_ll + v2_rr)

    p_avg = 0.5f0 * (p_ll + p_rr)
    T_ll = p_ll/(equations.R * rho_ll)
    T_rr = p_rr/(equations.R * rho_rr)
    T_avg = 0.5f0 * (T_ll + T_rr) * equations.c_v
    rho_e_ll = u_ll[4]
    rho_e_rr = u_rr[4]
    e_ll = rho_e_ll/rho_ll
    e_rr = rho_e_rr/rho_rr
    p_avg = 0.5f0 * (p_ll + p_rr)
    e_avg = 0.5f0 * (e_ll + e_rr)
    # Calculate fluxes depending on normal_direction
    f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
    f2 = f1 * v1_avg + p_avg * normal_direction[1] 
    f3 = f1 * v2_avg + p_avg * normal_direction[2]
    f4 = f1 * e_avg + p_avg * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) 
    return SVector(f1, f2, f3, f4, zero(eltype(u_ll)))
end



@inline function flux_lmars(u_ll, u_rr, normal_direction::AbstractVector,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
	a = 340
    # Unpack left and right state
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
    v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

    norm_ = norm(normal_direction)

    rho = 0.5f0 * (rho_ll + rho_rr)

    p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_rr - v_ll) / norm_
    v_interface = 0.5f0 * (v_ll + v_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

    if (v_interface > 0)
        f1, f2, f3, f4 = u_ll * v_interface
	f4 = f4 + p_ll * v_interface
    else
        f1, f2, f3, f4 = u_rr * v_interface
	f4 = f4 + p_rr * v_interface
    end

    return SVector(f1,
                   f2 + p_interface * normal_direction[1],
                   f3 + p_interface * normal_direction[2],
		   f4, zero(eltype(u_ll)))
end

@inline function Trixi.prim2cons(prim,
                                 equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, v1, v2, p, phi = prim
    rho_v1 = rho * v1
    rho_v2 = rho * v2
    rho_e = p/(equations.gamma-1)
    return SVector(rho, rho_v1, rho_v2, rho_e, phi)
end

@inline function cons2prim(u,
                           equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e, phi = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    p = (equations.gamma-1) * rho_e
    return SVector(rho, v1, v2, p, phi)
end

@inline function cons2cons(u, equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    return u
end

@inline function cons2entropy(u,
                              equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, rho_v1, rho_v2, rho_e = u

    w1 = log(rho_e*(equations.gamma-1)/rho^equations.gamma) - equations.gamma
     w4 = rho / (rho_e * (equations.gamma-1))

    return SVector(w1, zero(eltype(u)), zero(eltype(u)), w4, zero(eltype(u)))
end

@inline function entropy(cons,
                         equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    p = (equations.gamma-1) * cons[4]
    # Thermodynamic entropy
    s = log(p) - equations.gamma * log(cons[1])
    S = -s * cons[1] / (equations.gamma - 1)
    return S
end

@inline function pressure(cons,
                          equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    p = (equations.gamma-1) * cons[4]
    return p
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
                                     equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    elseif orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
                                     equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speed
    # left
    v_ll = (v1_ll * normal_direction[1]
	    + v2_ll * normal_direction[2])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
	    + v2_rr * normal_direction[2])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, orientation::Integer,
                               equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Get the velocity value in the appropriate direction
    if orientation == 1
        v_ll = v1_ll
        v_rr = v1_rr
    elseif orientation == 2
        v_ll = v2_ll
        v_rr = v2_rr
    end
    # Calculate sound speeds
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    return max(abs(v_ll) + c_ll, abs(v_rr) + c_rr)
end

# Less "cautious", i.e., less overestimating `λ_max` compared to `max_abs_speed_naive`
@inline function max_abs_speed(u_ll, u_rr, normal_direction::AbstractVector,
                               equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
    rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

    # Calculate normal velocities and sound speeds
    # left
    v_ll = (v1_ll * normal_direction[1]
	    + v2_ll * normal_direction[2])
    c_ll = sqrt(equations.gamma * p_ll / rho_ll)
    # right
    v_rr = (v1_rr * normal_direction[1]
	    + v2_rr * normal_direction[2])
    c_rr = sqrt(equations.gamma * p_rr / rho_rr)

    norm_ = norm(normal_direction)
    return max(abs(v_ll) + c_ll * norm_, abs(v_rr) + c_rr * norm_)
end

@inline function Trixi.max_abs_speeds(u,
                                equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
    rho, v1, v2, p = cons2prim(u, equations)
    c = sqrt(equations.gamma * p / rho)

    return abs(v1) + c, abs(v2) + c
end
end # @muladd
