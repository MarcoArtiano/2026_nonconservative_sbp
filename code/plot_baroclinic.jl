using CairoMakie
using LaTeXStrings
using DelimitedFiles

function plot_entropy(; basepath=pwd())

    filename = joinpath(basepath,
        "results/baroclinic",
        "data_5_16x8_20.0.dat")

    data = readdlm(filename, skipstart=1)


    fig = Figure(size=(1000, 500))
    kwargs = (xlabel=L"t \ \mathrm{[days]}", xlabelsize=30, ylabelsize=25, limits=((0, 20), nothing), xticklabelsize=17.0, yticklabelsize=17.0)
    ylabel = L"\underline{1} \underline{\underline{M}} \underline{U} - \underline{1} \underline{\underline{M}}\underline{U}_0"
    ylabel = L"\int_M U -\int_M U_0"
    ax1 = Axis(fig[1, 1]; ylabel=ylabel, xticks=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20], 
        kwargs...)

    lines!(ax1, data[:, 2] / 86400, data[:, end] .- data[1, end]; linewidth=3.5)
    save("results/baroclinic/entropy.pdf", fig)

end

plot_entropy()