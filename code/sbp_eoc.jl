using PrettyTables
using DelimitedFiles
using Trixi
using OrdinaryDiffEqSSPRK
using OrdinaryDiffEqLowStorageRK
using SummationByPartsOperators
using MuladdMacro

include("equations/sbp_volume.jl")
include("equations/compressible_euler_internal_energy_2d.jl")

function convergence_tests_curvilinear_2d_euler(; latex=false)

    accuracy_order = 4

    mesh_file_name = "mesh_convergence_test_deg2.mesh"

    num_nodes = [22, 44, 66, 88]
    @info "2D compressible Euler equations" accuracy_order
    _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes,
        accuracy_order, latex)

    accuracy_order = 6
    @info "2D compressible Euler equations" accuracy_order
    _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes,
        accuracy_order, latex)

    accuracy_order = 8
    @info "2D compressible Euler equations" accuracy_order
    _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes,
        accuracy_order, latex)

    return nothing
end

function _convergence_tests_curvilinear_2d_euler(; mesh_file_name, num_nodes, accuracy_order, latex=false)
    num_elements = Vector{Int}()
    errors = Vector{Float64}()

    for nnodes in num_nodes
        nelements = 16 # both meshes have fixed 16 elements so this value is just hard coded
        tol = 1.0e-13
        res = compute_errors_curvilinear_2d_euler(; mesh_file_name, nnodes,
            accuracy_order, tol)
        push!(num_elements, nelements)
        push!(errors, first(res.l2))
    end

    eoc = compute_eoc(num_nodes, errors)

    data = hcat(num_elements, num_nodes, errors, eoc)

    fmt = [
        fmt__printf("%3d", [1, 2]),
        fmt__printf("%.2e", [3]),
        fmt__printf("%.2f", [4])
    ]

    pretty_table(data)
    @show pretty_table(data)

    return nothing
end

function compute_errors_curvilinear_2d_euler(; mesh_file_name, nnodes,
    accuracy_order, tol)

    function initial_condition_convergence_shifted(x, t, equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
        c = 2
        A = 0.1
        L = sqrt(2)
        f = 1 / L
        ω = 2 * pi * f
        ini = c + A * sin(ω * (x[1] + x[2] - t))

        rho = ini
        rho_v1 = ini
        rho_v2 = ini
        rho_e = ini^2
        rho_internal = rho_e - 0.5f0 * rho_v1^2 / rho - 0.5f0 * rho_v2^2 / rho

        return SVector(rho, rho_v1, rho_v2, rho_internal, zero(eltype(x)))
    end
    @inline function source_terms_convergence_shifted(u, x, t,
        equations::CompressibleEulerInternalEnergyEquationsWithGravity2D)
        # Same settings as in `initial_condition`
        c = 2
        A = 0.1
        L = sqrt(2)
        f = 1 / L
        ω = 2 * pi * f
        γ = equations.gamma

        x1, x2 = x
        si, co = sincos(ω * (x1 + x2 - t))
        rho = c + A * si
        rho_x = ω * A * co
        # Note that d/dt rho = -d/dx rho = -d/dy rho.

        tmp = (2 * rho - 1) * (γ - 1)
        v1 = 1
        v2 = 1
        du1 = rho_x
        du2 = rho_x * (1 + tmp)
        du3 = du2
        du4 = 2 * rho_x * (rho + tmp)
        du_rhoe = du4 + (v1^2 + v2^2) / 2 * du1 - v1 * du2 - v2 * du3
        return SVector(du1, du2, du3, du_rhoe, zero(eltype(u)))
    end

    equations = CompressibleEulerInternalEnergyEquationsWithGravity2D(c_p=1004, c_v=717, gravity=9.81)
    initial_condition = initial_condition_convergence_shifted
    source_terms = source_terms_convergence_shifted

    D_SBP = derivative_operator(SummationByPartsOperators.MattssonAlmquistVanDerWeide2018Accurate(),
        derivative_order=1, accuracy_order=accuracy_order,
        xmin=-1.0, xmax=1.0, N=nnodes)

    solver = FDSBP(D_SBP, surface_integral=SurfaceIntegralStrongForm(flux_volume_ec_combined), volume_integral=VolumeIntegralFluxDifferencing(flux_volume_ec_combined))

    mesh_file = joinpath(@__DIR__, mesh_file_name)

    mesh = UnstructuredMesh2D(mesh_file, periodicity=true)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition,
        solver; source_terms)

    ode = semidiscretize(semi, (0.0, 2.0))

    summary_callback = SummaryCallback()

    analysis_interval = 100000
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    save_solution = SaveSolutionCallback(interval=10000,
        save_initial_solution=true,
        save_final_solution=true)
	stepsize_callback = StepsizeCallback(cfl = 0.01)
    sol = solve(ode,
        SSPRK43(thread=Trixi.True());
	    dt = 0.01,
        ode_default_options()...,
        abstol=tol, reltol=tol, callback=CallbackSet(stepsize_callback), adaptive = false)

    analysis_callback = AnalysisCallback(semi)
    return analysis_callback(sol)
end

function compute_eoc(Ns, errors)
    eoc = similar(errors)
    eoc[begin] = NaN # no EOC defined for the first grid
    for idx in Iterators.drop(eachindex(errors, Ns, eoc), 1)
        eoc[idx] = -(log(errors[idx] / errors[idx-1]) / log(Ns[idx] / Ns[idx-1]))
    end
    return eoc
end

convergence_tests_curvilinear_2d_euler()
