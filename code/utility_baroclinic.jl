# Unperturbed balanced steady-state.
# Returns primitive variables with only the velocity in longitudinal direction (rho, u, p).
# The other velocity components are zero.
function basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)
  # Parameters from Table 1 in the paper
  # Corresponding names in the paper are commented

  RealT = eltype(z)
  radius_earth = RealT(6.371229e6) # a
  half_width_parameter = 2f0           # b
  gravitational_acceleration = 9.81f0     # g
  k = 3f0           # k
  surface_pressure = RealT(1e5)         # p₀
  gas_constant = 287f0         # R
  surface_equatorial_temperature = 310.0f0       # T₀ᴱ
  surface_polar_temperature = 240.0f0       # T₀ᴾ
  lapse_rate = 0.005f0       # Γ
  angular_velocity = RealT(7.29212e-5)  # Ω

  # Distance to the center of the Earth
  r = z + radius_earth

  # In the paper: T₀
  temperature0 = 0.5f0 * (surface_equatorial_temperature + surface_polar_temperature)
  # In the paper: A, B, C, H
  const_a = 1 / lapse_rate
  const_b = (temperature0 - surface_polar_temperature) /
            (temperature0 * surface_polar_temperature)
  const_c = 0.5f0 * (k + 2) * (surface_equatorial_temperature - surface_polar_temperature) /
            (surface_equatorial_temperature * surface_polar_temperature)
  const_h = gas_constant * temperature0 / gravitational_acceleration

  # In the paper: (r - a) / bH
  scaled_z = z / (half_width_parameter * const_h)

  # Temporary variables
  temp1 = exp(lapse_rate / temperature0 * z)
  temp2 = exp(-scaled_z^2)

  # In the paper: ̃τ₁, ̃τ₂
  tau1 = const_a * lapse_rate / temperature0 * temp1 +
         const_b * (1 - 2 * scaled_z^2) * temp2
  tau2 = const_c * (1 - 2 * scaled_z^2) * temp2

  # In the paper: ∫τ₁(r') dr', ∫τ₂(r') dr'
  inttau1 = const_a * (temp1 - 1) + const_b * z * temp2
  inttau2 = const_c * z * temp2

  # Temporary variables
  temp3 = r / radius_earth * cos(lat)
  temp4 = temp3^k - k / (k + 2) * temp3^(k + 2)

  # In the paper: T
  temperature = 1 / ((r / radius_earth)^2 * (tau1 - tau2 * temp4))

  # In the paper: U, u (zonal wind, first component of spherical velocity)
  big_u = gravitational_acceleration / radius_earth * k * temperature * inttau2 *
          (temp3^(k - 1) - temp3^(k + 1))
  temp5 = radius_earth * cos(lat)
  u = -angular_velocity * temp5 + sqrt(angular_velocity^2 * temp5^2 + temp5 * big_u)

  # Hydrostatic pressure
  p = surface_pressure *
      exp(-gravitational_acceleration / gas_constant * (inttau1 - inttau2 * temp4))

  # Density (via ideal gas law)
  rho = p / (gas_constant * temperature)

  return rho, u, p
end

# Perturbation as in Equations 25 and 26 of the paper (analytical derivative)
function perturbation_stream_function(lon, lat, z)
  # Parameters from Table 1 in the paper
  # Corresponding names in the paper are commented
  perturbation_radius = 1 / 6f0      # d₀ / a
  perturbed_wind_amplitude = 1f0     # Vₚ
  perturbation_lon = pi / 9f0     # Longitude of perturbation location
  perturbation_lat = 2 * pi / 9f0 # Latitude of perturbation location
  pertz = 15000f0    # Perturbation height cap

  # Great circle distance (d in the paper) divided by a (radius of the Earth)
  # because we never actually need d without dividing by a
  great_circle_distance_by_a = acos(sin(perturbation_lat) * sin(lat) +
                                    cos(perturbation_lat) * cos(lat) *
                                    cos(lon - perturbation_lon))

  # In the first case, the vertical taper function is per definition zero.
  # In the second case, the stream function is per definition zero.
  if z > pertz || great_circle_distance_by_a > perturbation_radius
    return 0, 0
  end

  # Vertical tapering of stream function
  perttaper = 1 - 3f0 * z^2 / pertz^2 + 2 * z^3 / pertz^3

  # sin/cos(pi * d / (2 * d_0)) in the paper
  sin_, cos_ = sincos(0.5f0 * pi * great_circle_distance_by_a / perturbation_radius)

  # Common factor for both u and v
  factor = 16f0 / (3f0 * sqrt(3f0)) * perturbed_wind_amplitude * perttaper * cos_^3 * sin_

  u_perturbation = -factor * (-sin(perturbation_lat) * cos(lat) +
                              cos(perturbation_lat) * sin(lat) * cos(lon - perturbation_lon)) /
                   sin(great_circle_distance_by_a)

  v_perturbation = factor * cos(perturbation_lat) * sin(lon - perturbation_lon) /
                   sin(great_circle_distance_by_a)

  return u_perturbation, v_perturbation
end

function cartesian_to_sphere(x)
  r = Trixi.norm(x)
  lambda = atan(x[2], x[1])
  if lambda < 0
    lambda += 2 * pi
  end
  phi = asin(x[3] / r)

  return lambda, phi, r
end

using PyPlot

matplotlib = PyPlot.matplotlib
rcParams = matplotlib["rcParams"]

rcParams["pdf.fonttype"] = 42  # Usa Type 42 fonts TrueType, che vengono embedded subset nel PDF
rcParams["ps.fonttype"] = 42

rcParams["text.usetex"] = true  # Se vuoi usare LaTeX per i testi, altrimenti false
rcParams["font.family"] = "serif"  #

@inline function contour_baroclinic(sol, semi, (Kh, Kv), nvisnodes, equations_euler, Tf, surface_flux)

  data_layer_euler, nodes_layer_euler = retrieve_values_per_layer(Kh, Kv, semi, sol)

  p_euler = retrieve_pressure_layer(data_layer_euler, equations_euler)
  T_euler = retrieve_temperature_layer(data_layer_euler, equations_euler)

  p_euler, x, y, z = spectral_interpolation3D(p_euler, nodes_layer_euler, semi, nvisnodes)

  T_euler, _, _, _ = spectral_interpolation3D(T_euler, nodes_layer_euler, semi, nvisnodes)

  lon, lat = cart2sphere(x, y, z)

  lonshift = 60
  lon = @. mod(lon + 180 - lonshift, 360) - 180

  mask = (lat .>= 0) .& (lon .>= -lonshift)

  lon = lon[mask]
  lat = lat[mask]
  p_euler = p_euler[mask]
  t_euler = T_euler[mask]
  plevels = vcat([955], [960 + 5i for i ∈ 0:11], [1020])
  plevels = 11
  pnorm = matplotlib.colors.TwoSlopeNorm(vmin=955, vcenter=990, vmax=1025)
  pnorm = nothing
  cmap = ColorMap("plasma").copy()
  shrinkcb = 0.7
  fig, axs = subplots(1, 1, figsize=(27, 20))
  cs1 = axs.tricontour(lon, lat, p_euler; levels=plevels, colors=("k",))
  cset = axs.tricontourf(
    lon,
    lat,
    p_euler;
    levels=plevels,
    cmap,
    norm=pnorm,
    extend="neither",
  )

  axs[:set_title](raw"Day 10, Surface Pressure [hPa]", fontsize=40)

  cbar_ax = fig.add_axes([0.94, 0.25, 0.015, 0.5])  # Cambia qui per regolare posizione
  cbar = colorbar(
    cset,
    cax=cbar_ax,
    ax=axs,
    ticks=plevels isa Int ? nothing : plevels[1+2:2:end-1],
    shrink=shrinkcb)
  plt.subplots_adjust(wspace=0.05)
  axs.tick_params(axis="both", labelsize=40)  # Sostituisci 14 con la dimensione desiderata

  xticks = [-60, -30, 0, 30, 60, 90, 120, 150, 180]
  xticklabels =
    ["0", "30E", "60E", "90E", "120E", "150E", "180", "150W", "120W"]
  yticks = [0, 30, 60, 90]
  yticklabels = ["0", "30N", "60N", "90N"]
  axs.set_xticks(xticks)
  axs.set_xticklabels(xticklabels)
  axs.set_yticks(yticks)
  axs.set_yticklabels(yticklabels)
  axs.set_xlim([xticks[1], xticks[end]])
  axs.set_ylim([yticks[1], yticks[end]])
  axs.set_aspect(1)
  cbar[:ax][:tick_params](labelsize=30)

  PyPlot.savefig(joinpath("results/baroclinic/contour_pressure_euler_$(Kh)x$(Kv)_$(Tf)_$(surface_flux).png"), dpi=300, bbox_inches="tight")
  close(fig)

  fig, axs = subplots(1, 1, figsize=(27, 20))
  cs1 = axs.tricontour(lon, lat, t_euler; levels=plevels, colors=("k",))
  cset = axs.tricontourf(
    lon,
    lat,
    t_euler;
    levels=plevels,
    cmap,
    norm=pnorm,
    extend="neither",
  )

  axs[:set_title](raw"Day 10, Temperature [K]", fontsize=40)

  cbar_ax = fig.add_axes([0.94, 0.25, 0.015, 0.5])  # Cambia qui per regolare posizione
  cbar = colorbar(
    cset,
    cax=cbar_ax,
    ax=axs,
    ticks=plevels isa Int ? nothing : plevels[1+2:2:end-1],
    shrink=shrinkcb)
  plt.subplots_adjust(wspace=0.05)
  axs.tick_params(axis="both", labelsize=40)  # Sostituisci 14 con la dimensione desiderata

  xticks = [-60, -30, 0, 30, 60, 90, 120, 150, 180]
  xticklabels =
    ["0", "30E", "60E", "90E", "120E", "150E", "180", "150W", "120W"]
  yticks = [0, 30, 60, 90]
  yticklabels = ["0", "30N", "60N", "90N"]
  axs.set_xticks(xticks)
  axs.set_xticklabels(xticklabels)
  axs.set_yticks(yticks)
  axs.set_yticklabels(yticklabels)
  axs.set_xlim([xticks[1], xticks[end]])
  axs.set_ylim([yticks[1], yticks[end]])
  axs.set_aspect(1)
  cbar[:ax][:tick_params](labelsize=30)

  PyPlot.savefig(joinpath("results/baroclinic/contour_temperature_euler_$(Kh)x$(Kv)_$(Tf)_$(surface_flux).png"), dpi=300, bbox_inches="tight")
  close(fig)


  fig, axs = subplots(1, 1, figsize=(27, 20))
  cs1 = axs.tricontour(lon, lat, p_euler; levels=plevels, colors=("k",))
  cset = axs.tricontourf(
    lon,
    lat,
    p_euler;
    levels=plevels,
    cmap,
    norm=pnorm,
    extend="neither",
  )

end

function retrieve_pressure_layer(data_layer, equations)
  size_ = size(data_layer)
  p = zeros(Float64, size_[2:end]...)
  for i in 1:size_[2]
    for j in 1:size_[3]
      for element in 1:size_[4]
        p[i, j, element] = cons2pressure(data_layer[:, i, j, element], equations)
      end
    end
  end
  return p

end

function retrieve_temperature_layer(data_layer, equations)
  size_ = size(data_layer)
  T = zeros(Float64, size_[2:end]...)
  for i in 1:size_[2]
    for j in 1:size_[3]
      for element in 1:size_[4]
        T[i, j, element] = cons2temperature(data_layer[:, i, j, element], equations)
      end
    end
  end
  return T

end

function cons2temperature(u, equations::CompressibleEulerInternalEnergyEquationsWithGravity3D)
  rho, rho_v1, rho_v2, rho_v3, rho_e = u

  p = (equations.gamma - 1) * rho_e
  return p / (rho * 287)
end

function cons2pressure(u, equations::CompressibleEulerInternalEnergyEquationsWithGravity3D)
  rho, rho_v1, rho_v2, rho_v3, rho_e = u

  p = (equations.gamma - 1) * rho_e
  return p / 1e2
end


function cons2temperature(u, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
  rho, rho_v1, rho_v2, rho_v3, rho_theta = u

  p = equations.K * rho_theta^equations.gamma
  return p / (rho * 287)
end

function cons2pressure(u, equations::CompressibleEulerPotentialTemperatureEquationsWithGravity3D)
  rho, rho_v1, rho_v2, rho_v3, rho_theta = u

  p = equations.K * rho_theta^equations.gamma
  return p / 1e2
end

function cart2sphere(x, y, z)
  r = similar(x)
  r .= sqrt.(x .^ 2 + y .^ 2 + z .^ 2)

  lat = similar(x)
  lon = similar(x)

  lat = asin.(z ./ r) .* (180 / π)
  lon = atan.(y, x) .* (180 / π)

  return vec(lon), vec(lat)
end

@inline function retrieve_values_per_layer(Kh, Kv, semi, sol; layer=1)

  node_coordinates = semi.cache.elements.node_coordinates
  data = Trixi.wrap_array(sol.u[end], semi)
  size_sol = size(data)
  size_coords = size(node_coordinates)
  sphere_layer = (Kh^2*(layer-1)+1):Kh^2*Kv:6*Kv*Kh^2
  data_layer = zeros(Float64, size_sol[1:end-2]..., Kh^2 * 6)
  nodes_layer = zeros(Float64, size_coords[1:end-2]..., Kh^2 * 6)

  for block in 1:6
    data_layer[:, :, :, (Kh^2*(block-1)+1):Kh^2*block] .= data[:, :, :, 1, sphere_layer[block]:(sphere_layer[block]+Kh^2-1)]
    nodes_layer[:, :, :, (Kh^2*(block-1)+1):Kh^2*block] .= node_coordinates[:, :, :, 1, sphere_layer[block]:(sphere_layer[block]+Kh^2-1)]
  end
  return data_layer, nodes_layer
end

function spectral_interpolation3D(data_layer, nodes_layer, semi, nvisnodes)
  nvars = 1
  size_ = size(nodes_layer)
  n_nodes_2d = Trixi.nnodes(semi.solver)^2
  n_elements = size_[end]
  #  plotting_interp_matrix = plotting_interpolation_matrix_no_boundary(semi.solver; nvisnodes=nvisnodes)
  plotting_interp_matrix = Trixi.plotting_interpolation_matrix(semi.solver; nvisnodes=nvisnodes)

  uEltype = eltype(data_layer)
  x = reshape(view(nodes_layer, 1, :, :, :), n_nodes_2d,
    n_elements)
  y = reshape(view(nodes_layer, 2, :, :, :), n_nodes_2d,
    n_elements)
  z = reshape(view(nodes_layer, 3, :, :, :), n_nodes_2d,
    n_elements)

  u_extracted = StructArray{SVector{nvars,uEltype}}(ntuple(_ -> similar(x,
      (n_nodes_2d,
        n_elements)),
    nvars))
  for element in 1:n_elements
    sk = 1
    for j in eachnode(semi.solver), i in eachnode(semi.solver)
      u_node = SVector(data_layer[i, j, element])
      u_extracted[sk, element] = u_node
      sk += 1
    end
  end
  uplot = StructArray{SVector{nvars,uEltype}}(map(x -> plotting_interp_matrix * x,
    StructArrays.components(u_extracted)))
  xplot, yplot, zplot = plotting_interp_matrix * x, plotting_interp_matrix * y, plotting_interp_matrix * z

  uplot = StructArray{SVector{nvars,uEltype}}(map(x -> plotting_interp_matrix * x,
    StructArrays.components(u_extracted)))

  return getindex.(vec(uplot), 1), vec(xplot), vec(yplot), vec(zplot)
end

function plotting_interpolation_matrix_no_boundary(dg::DGSEM;
  nvisnodes=2 * length(dg.basis.nodes))

  dξ = 2 / nvisnodes
  interia = [-1 + (j - 1 / 2) * dξ for j ∈ 1:nvisnodes]
  Vp1D = Trixi.polynomial_interpolation_matrix(dg.basis.nodes, interia)
  # For quadrilateral elements, interpolation to plotting nodes involves applying a 1D interpolation
  # operator to each line of nodes. This is equivalent to multiplying the vector containing all node
  # node coordinates on an element by a Kronecker product of the 1D interpolation operator (e.g., a
  # multi-dimensional interpolation operator).
  return Trixi.kron(Vp1D, Vp1D)
end
