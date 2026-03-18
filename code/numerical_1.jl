using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqSSPRK
using Trixi
using Trixi: @muladd
using Plots

include("equations/utilities.jl")
include("equations/advection_variable_coefficient.jl")
equations = AdvectionVariableCoefficientEquation1D()

function initial_condition(x, t, equation::AdvectionVariableCoefficientEquation1D)
    RealT = eltype(x)
    k = 1
    a = 2 + cospi(x[1])
    return SVector(2 + sinpi(k * (x[1] - convert(RealT, 0.7))), a)
end

function run_1(polydeg)
    volume_flux = (flux_zero, flux_entropy_conservative)
    solver = DGSEM(polydeg=polydeg, surface_flux=(flux_zero, flux_entropy_conservative),
        volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

    coordinates_min = -1.0
    coordinates_max = 1.0
    mesh = TreeMesh(coordinates_min, coordinates_max,
        initial_refinement_level=5,
        n_cells_max=10_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition,
        solver)

    ###############################################################################
    # ODE solvers, callbacks etc.
    tspan = (0.0, 1.0)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 100
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
        extra_analysis_errors=(:l2_error_primitive,
                :linf_error_primitive), extra_analysis_integrals = (entropy, ), save_analysis = true,
                output_directory = pwd() * "/results/numerical_1/",
                analysis_filename = "data_$(polydeg).dat",)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)
        stepsize_callback = StepsizeCallback(cfl = 0.01)
    callbacks = CallbackSet(summary_callback,
        analysis_callback, alive_callback, stepsize_callback)


    ###############################################################################
    # run the simulation

    sol = solve(ode,
	    SSPRK43(thread=Trixi.True());
	    dt = 0.1,
        ode_default_options()..., callback=callbacks, adaptive = false);
    return sol
end

polydegs = (1,2,3,4,5)
for polydeg in polydegs
	run_1(polydeg)
end

@info "Detailed raw data saved in results/numerical_1/data_*.dat"

using DelimitedFiles

function compute_table(polydegs; basepath=pwd())
    table = zeros(length(polydegs), 3)

    for (i, p) in enumerate(polydegs)
        filename = joinpath(basepath, "results/numerical_1", "data_$(p).dat")
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

@info "Results: polynomial degree, fully-discrete entropy change, semi-disrete entropy change"
table = compute_table(polydegs)'
