using OrdinaryDiffEqLowStorageRK
using OrdinaryDiffEqSSPRK
using Trixi
using Trixi: @muladd
using Plots
using Optim

include("equations/polynomial_equations.jl")
###############################################################################

function initial_condition_linear_stability(x, t, equation::PolynomialEquation1D)
    RealT = eltype(x)
    #k = 1
    #return SVector(2 + sinpi(k * (x[1] - convert(RealT, 0.7))))
    return SVector(sinpi(pi * x[1]))
end

function run_2(polydeg, m, n, bool)

    alpha = (m + 1) / (m + n + 1)
    equations = PolynomialEquation1D(m, n, alpha)
    if bool == true
        volume_flux = (flux_zero, flux_entropy_general)
        solver = DGSEM(polydeg=polydeg, surface_flux=(flux_zero, flux_entropy_general),
            volume_integral=VolumeIntegralFluxDifferencing(volume_flux))
        name = "noncons"
    else
        volume_flux = (flux_conservative_even, flux_nonconservative_even)
        solver = DGSEM(polydeg=polydeg, surface_flux=(flux_conservative_even, flux_nonconservative_even),
            volume_integral=VolumeIntegralFluxDifferencing(volume_flux))
        name = "cons"
    end
    coordinates_min = -1.0
    coordinates_max = 1.0
    mesh = TreeMesh(coordinates_min, coordinates_max,
        initial_refinement_level=5,
        n_cells_max=10_000)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_linear_stability,
        solver)
    r = m + n
    f(x) = n * (r - 1) * pi * sinpi(x)^(r - 2) * cospi(x)
    res = optimize(f, -1.0, 1.0)
    fmin = Optim.minimum(res)
    T_max = -1 / fmin
    T = T_max * 0.5
    ###############################################################################
    # ODE solvers, callbacks etc.
    tspan = (0.0, T)
    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 100
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
        extra_analysis_errors=(:l2_error_primitive,
            :linf_error_primitive), extra_analysis_integrals=(entropy, velocity), save_analysis=true,
        output_directory=pwd() * "/results/numerical_2/",
        analysis_filename="data_$(polydeg)_$(m)_$(n)_$(bool).dat",)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    stepsize_callback = StepsizeCallback(cfl=0.001)

    callbacks = CallbackSet(summary_callback,
        analysis_callback, alive_callback,
        stepsize_callback
    )

    ###############################################################################
    # run the simulation

    sol = solve(ode, SSPRK43();
        dt=0.01, # solve needs some value here but it will be overwritten by the stepsize_callback
        ode_default_options()..., callback=callbacks, adaptive=false)

    a = plot(sol)
end

ms = (2, 3, 4, 5)
ns = (2, 3, 4, 5)
polydeg = 2
polydegs = (1, 2, 3, 4, 5)
bools = (true, false)
for polydeg in polydegs
    for m in ms
        for n in ns
            for vbool in bools
                run_2(polydeg, m, n, vbool)
            end
        end
    end
end


using Printf

using DelimitedFiles

function compute_tables(polydegs, ms, ns, bools; basepath=pwd())

    tables = Dict()

    for m in ms, n in ns, vbool in bools

        key = (m, n, vbool)

        table = zeros(length(polydegs), 4)
        row = 1

        for p in polydegs

            filename = joinpath(basepath,
                "results/numerical_2",
                "data_$(p)_$(m)_$(n)_$(vbool).dat")

            data = readdlm(filename, skipstart=1)

            U = data[:, end-1]
            wMdt = data[:, end-2]
            u = data[:, end]
            Umean = abs(U[end] - U[1]) / U[1]
            umean = u[end]

            table[row, :] .= (p, Umean, wMdt[end], umean)
            row += 1
        end

        tables[key] = table'
    end

    return tables
end

table = compute_tables(polydegs, ms, ns, bools)
