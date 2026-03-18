using Trixi
using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqSSPRK
using MuladdMacro
using TrixiAtmo
using DelimitedFiles
#using LaTeX
###############################################################################
# Setup for the baroclinic instability test
include("equations/compressible_euler_internal_energy_3d.jl")
include("utility_baroclinic.jl")
# Initial condition for an idealized baroclinic instability test
# https://doi.org/10.1002/qj.2241, Section 3.2 and Appendix A
function initial_condition_baroclinic_instability(x, t,
    equations::CompressibleEulerInternalEnergyEquationsWithGravity3D)
    lon, lat, r = cartesian_to_sphere(x)
    RealT = eltype(x)
    radius_earth = RealT(6.371229e6)
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Stream function type perturbation
    u_perturbation, v_perturbation = perturbation_stream_function(lon, lat, z)

    u += u_perturbation
    v = v_perturbation

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u - sin(lat) * cos(lon) * v
    v2 = cos(lon) * u - sin(lat) * sin(lon) * v
    v3 = cos(lat) * v
    radius_earth = RealT(6.371229e6)  # a
    gravitational_acceleration = RealT(9.81)    # g

    r = Trixi.norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0)
    if z > 0
        r = Trixi.norm(x)
    else
        r = -(2 * radius_earth^3) / (x[1]^2 + x[2]^2 + x[3]^2)
    end
    r = -Trixi.norm(x)
    phi = radius_earth^2 * gravitational_acceleration / r

    return prim2cons(SVector(rho, v1, v2, v3, p, phi), equations)
end

# Steady state for RHS correction below
function steady_state_baroclinic_instability(x, t, equations::CompressibleEulerInternalEnergyEquationsWithGravity3D)
    lon, lat, r = cartesian_to_sphere(x)
    RealT = eltype(x)
    radius_earth = RealT(6.371229e6)
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u
    v2 = cos(lon) * u
    v3 = 0
    radius_earth = RealT(6.371229e6)  # a
    gravitational_acceleration = RealT(9.81)     # g

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0)

    if z > 0
        r = norm(x)
    else
        r = -(2 * radius_earth^3) / (x[1]^2 + x[2]^2 + x[3]^2)
    end
    r = -norm(x)
    phi = radius_earth^2 * gravitational_acceleration / r

    return prim2cons(SVector(rho, v1, v2, v3, p, phi), equations)
end


@inline function source_terms_baroclinic_instability(u, x, t,
    equations::CompressibleEulerInternalEnergyEquationsWithGravity3D)

    RealT = eltype(u)

    radius_earth = RealT(6.371229e6)  # a
    gravitational_acceleration = RealT(9.81)     # g
    angular_velocity = RealT(7.29212e-5)  # Ω

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0)
    r = z + radius_earth

    du1 = zero(eltype(u))
    du2 = zero(eltype(u))
    du3 = zero(eltype(u))
    du4 = zero(eltype(u))
    du5 = zero(eltype(u))
    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 -= -2 * angular_velocity * u[3]
    du3 -= 2 * angular_velocity * u[2]

    return SVector(du1, du2, du3, du4, du5, zero(eltype(u)))
end

###############################################################################
# Start of the actual elixir, semidiscretization of the problem

function main(equations, surface_flux, volume_flux, T, filename, trees_per_cube_face, polydeg, RealT)

    initial_condition = initial_condition_baroclinic_instability

    boundary_conditions = Dict(:inside => boundary_condition_slip_wall,
        :outside => boundary_condition_slip_wall)

    solver = DGSEM(RealT=RealT, polydeg=polydeg, surface_flux=surface_flux,
        volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

    mesh = Trixi.P4estMeshCubedSphere(trees_per_cube_face..., RealT(6.371229e6), RealT(30000.0), RealT=RealT,
        polydeg=polydeg, initial_refinement_level=0)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
        source_terms=source_terms_baroclinic_instability,
        boundary_conditions=boundary_conditions)

    ###############################################################################
    # ODE solvers, callbacks etc.
    tspan = (0.0, RealT(T * 24 * 60 * 60.0)) # time in seconds for 10 days
    #tspan = (0.0, 1000*dt) # time in seconds for 10 days

    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 100
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval, extra_analysis_integrals=(entropy,), save_analysis=true,
        output_directory=joinpath(@__DIR__, "results", "baroclinic"),
        analysis_filename="data_$(polydeg)_$(trees_per_cube_face[1])x$(trees_per_cube_face[2])_$(T).dat",)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    save_solution = SaveSolutionCallback(interval=400, save_initial_solution=true,
        save_final_solution=true,
        output_directory="results/baroclinic/" * filename * "$(trees_per_cube_face[1])x$(trees_per_cube_face[2])_$(polydeg)")

    callbacks = CallbackSet(summary_callback,
        analysis_callback,
        alive_callback,
        save_solution
    )

    tol = 1e-6
    ###############################################################################
    # Use a Runge-Kutta method with automatic (error based) time step size control
    # Enable threading of the RK method for better performance on multiple threads
    sol = solve(ode,
        RDPK3SpFSAL49(thread=Trixi.True());
        dt=RealT(1.0),
        abstol=tol,
        reltol=tol, ode_default_options()...,
        callback=callbacks, adaptive=true)
    return sol, semi
end

function run_baroclinic(; T, trees_per_cube_face=(8, 4), polydeg=5, RealT=Float64)
    surface_flux = flux_surface_combined
    volume_flux = flux_volume_ec_combined
    equations = CompressibleEulerInternalEnergyEquationsWithGravity3D(c_p=1004,
        c_v=717,
        gravity=9.81)
    equations_euler = equations
    filename = "data_" * string(T)
    sol, semi = main(equations, surface_flux, volume_flux, T, filename, trees_per_cube_face, polydeg, RealT)

    contour_baroclinic(sol, semi, trees_per_cube_face, 2 * (polydeg + 1), equations_euler, T, surface_flux)
end

run_baroclinic(T = 10.0, trees_per_cube_face = (16, 8))

run_baroclinic(T = 20.0, trees_per_cube_face = (16, 8))

