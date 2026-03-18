using MuladdMacro
using Trixi
using TrixiShallowWater

include("equations/escalante_saintemarie_1d.jl")

using OrdinaryDiffEqSSPRK

###############################################################################
# semidiscretization of the shallow water equations with a discontinuous
# bottom topography function for a fully wet configuration

equations = EscalanteSainteMarieEquations1D(gravity=1.0, b0=0.1)
alpha1 = 1/2
alpha2 = 1
alpha3 = 2/3

function initial_condition_periodic(x, t, equations::EscalanteSainteMarieEquations1D)
    h = 1 + exp(sinpi(2 * x[1]))
    v1 = 1
    v2 = 1
    p = 10
    b = 0.1*exp(sinpi(2 * x[1]))
    H = h+b
    return prim2cons(SVector(H, v1, v2, p, b), equations)
end

initial_condition = initial_condition_periodic
###############################################################################
# Get the DG approximation space
function run_1(polydeg)
volume_flux = (flux_conservative, flux_nonconservative)
surface_flux = (flux_conservative, flux_nonconservative)
solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux,
    volume_integral=VolumeIntegralFluxDifferencing(volume_flux))
boundary_condition = BoundaryConditionDirichlet(initial_condition)
###############################################################################
# Get the TreeMesh and setup a periodic mesh
coordinates_min = 0.0
coordinates_max = 1.0
mesh = TreeMesh(coordinates_min, coordinates_max,
    initial_refinement_level=7,
    n_cells_max=10_000,
    periodicity=true)

# Create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms = source_terms_escalante)

###############################################################################
# ODE solver
tspan = (0.0, 0.1)
ode = semidiscretize(semi, tspan)

###############################################################################
# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval, extra_analysis_integrals = (entropy, ), save_analysis = true,
			output_directory = pwd() * "/results/escalante_sainte_marie/",
			analysis_filename = "data_$(polydeg).dat",)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

stepsize_callback = StepsizeCallback(cfl=0.1)

callbacks = CallbackSet(
    summary_callback,
    analysis_callback,
    alive_callback,
    stepsize_callback
)

###############################################################################
# run the simulation

sol = solve(ode,
    SSPRK43();
    dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
    ode_default_options()..., callback=callbacks, adaptive=false);
end


polydegs = (1,2,3,4,5)
for polydeg in polydegs
	run_1(polydeg)
end

using DelimitedFiles

function compute_table(polydegs; basepath=pwd())
    table = zeros(length(polydegs), 3)

    for (i, p) in enumerate(polydegs)
        filename = joinpath(basepath, "results/escalante_sainte_marie", "data_$(p).dat")
        data = readdlm(filename, skipstart=1)

        U = data[:, end]
        wMdt = data[:, end-1]

        Umean = abs(U[end] - U[1]) / U[1]

        table[i, 1] = p
        table[i, 2] = Umean
        table[i, 3] = wMdt[end]   # oppure mean(wMdt)
    end

    return table
end

table = compute_table(polydegs)'
