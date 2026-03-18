using MuladdMacro
using Trixi
using TrixiShallowWater
using Plots

include("equations/escalante_saintemarie_2d.jl")

using OrdinaryDiffEqSSPRK
using Trixi

equations = EscalanteSainteMarieEquations2D(gravity=9.81, b0=0.1, H0 = 3.0)
alpha1 = 1/2
alpha2 = 1
alpha3 = 2/3

# An initial condition with constant total water height and zero velocities to test well-balancedness.
# Note, this routine is used to compute errors in the analysis callback but the initialization is
# overwritten by `initial_condition_discontinuous_well_balancedness` below.
function initial_condition_well_balancedness(x, t, equations::EscalanteSainteMarieEquations2D)
    # Calculate primitive variables
    H = equations.H0
    v1 = 0.0
    v2 = 0.0

    x1, x2 = x
    b = (1.5 / exp(0.5 * ((x1 - 1.0)^2 + (x2 - 1.0)^2)) +
         0.75 / exp(0.5 * ((x1 + 1.0)^2 + (x2 + 1.0)^2)))

    p = 0.0
    v3 = 0.0

    return prim2cons(SVector(H, v1, v2, v3, p, b), equations)
end

initial_condition = initial_condition_well_balancedness

boundary_condition = (; all = boundary_condition_slip_wall)

###############################################################################
# Get the DG approximation space

volume_flux = (flux_conservative, flux_nonconservative)
surface_flux = (flux_conservative, flux_nonconservative)
function run_wb(polydeg)
    # Create the solver
    solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux,
                volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
    mesh_file_name = "mesh_convergence_test_deg2.mesh"
    mesh_file = joinpath(@__DIR__, mesh_file_name)

    mesh = UnstructuredMesh2D(mesh_file, periodicity=true)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition,
            solver, source_terms = source_terms_escalante)

    ###############################################################################
    # ODE solvers, callbacks, etc.

    tspan = (0.0, 100.0)
    ode = semidiscretize(semi, tspan)
    function initial_condition_discontinuous_well_balancedness(x, t, element_id,
                                                            equations::EscalanteSainteMarieEquations2D)
        # Set the background values
        H = equations.H0
        v1 = 0.0
        v2 = 0.0
        v3 = 0.0
        b = 0.0
        p = 0.0
        # Setup a discontinuous bottom topography using the element id number
        if element_id == 7
            b = 2.0 + 0.5 * sin(2.0 * pi * x[1]) + 0.5 * cos(2.0 * pi * x[2])
        end

        return prim2cons(SVector(H, v1, v2, v3, p, b), equations)
    end

    # point to the data we want to augment
    u = Trixi.wrap_array(ode.u0, semi)
    # reset the initial condition
    for element in eachelement(semi.solver, semi.cache)
        for j in eachnode(semi.solver), i in eachnode(semi.solver)
            x_node = Trixi.get_node_coords(semi.cache.elements.node_coordinates, equations,
                                        semi.solver, i, j, element)
            u_node = initial_condition_discontinuous_well_balancedness(x_node, first(tspan),
                                                                    element, equations)
            Trixi.set_node_vars!(u, u_node, equations, semi.solver, i, j, element)
        end
    end
    ###############################################################################

    summary_callback = SummaryCallback()

    analysis_interval = 10000
    analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                        extra_analysis_errors = (:conservation_error,),
                                        extra_analysis_integrals = (lake_at_rest_error, velocity_1, velocity_2, velocity_3, pressure_), save_analysis = true, output_directory = joinpath(@__DIR__, "results", "escalante_sainte_marie"),
                analysis_filename = "data_wb_$(polydeg).dat",)

    alive_callback = AliveCallback(analysis_interval = analysis_interval)

    stepsize_callback = StepsizeCallback(cfl = 1.0)

    callbacks = CallbackSet(summary_callback,
                            analysis_callback,
                            alive_callback,
                            stepsize_callback)

    ###############################################################################
    # run the simulation
    # sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
    sol = solve(ode, SSPRK43(thread=Trixi.True()),
                dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
                adaptive = false,
                save_everystep = false, callback = callbacks, maxiters = 2e5);
    return sol
end

polydegs = (2,3,4,5)
for polydeg in polydegs
	run_wb(polydeg)
end

@info "Detailed raw data saved in results/escalante_sainte_marie/data_wb_*.dat"

using DelimitedFiles

function compute_table(polydegs; basepath=@__DIR__)
    table = zeros(length(polydegs), 5)

    for (i, p) in enumerate(polydegs)
        filename = joinpath(basepath, "results/escalante_sainte_marie", "data_wb_$(p).dat")
        data = readdlm(filename, skipstart=1)

        H = data[:, end-4]
        v1 = data[:, end-3]
        v2 = data[:, end-2]
        v3 = data[:, end-1]
        p = data[:,end]


        table[i, 1] = abs(H[end])
        table[i, 2] = abs(v1[end])
        table[i, 3] = abs(v2[end])
        table[i, 4] = abs(v3[end])
        table[i, 5] = abs(p[end])   # oppure mean(wMdt)
    end

    return table
end

@info "Results per polynomial degree: lake at rest error, velocity v1, velocity v2, velocity w, pressure p"
table = compute_table(polydegs)'