using Trixi: False, VolumeIntegralFluxDifferencing, FDSBP, eachelement, eachnode, get_contravariant_vector, @threaded, nnodes, multiply_add_to_node_vars!, AbstractDerivativeOperator, get_surface_normal, True, combine_conservative_and_nonconservative_fluxes, eachinterface, get_one_sided_surface_node_vars, get_node_coords


function Trixi.calc_volume_integral!(du, u,
    mesh::UnstructuredMesh2D,
    have_nonconservative_terms::False, equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::FDSBP, cache)

    D = cache.derivative_split # SBP derivative operator
    @unpack contravariant_vectors = cache.elements
    @unpack volume_flux = volume_integral

    if nvariables(equations) == 1
        u_vectors = reshape(reinterpret(SVector{nvariables(equations),eltype(u)}, u),
            nnodes(dg), nnodes(dg), nelements(dg, cache))
    else
        u_vectors = reinterpret(reshape, SVector{nvariables(equations),eltype(u)}, u)
    end

    @threaded for element in eachelement(dg, cache)
        u_element = view(u_vectors, :, :, element)

        @inbounds for j in eachnode(dg), i in eachnode(dg)
            u_node = u_element[i, j]
            Ja1_node = get_contravariant_vector(1, contravariant_vectors, i, j, element)
            Ja2_node = get_contravariant_vector(2, contravariant_vectors, i, j, element)

            for ii in (i+1):nnodes(dg)
                d_i_ii = D[i, ii]
                d_ii_i = D[ii, i]

                if iszero(d_i_ii) && iszero(d_ii_i)
                    continue
                end

                u_node_ii = u_element[ii, j]
                Ja1_node_ii = get_contravariant_vector(1, contravariant_vectors, ii, j, element)
                Ja1_avg = (Ja1_node + Ja1_node_ii) * 0.5f0

                fluxtilde1 = volume_flux(u_node, u_node_ii, Ja1_avg, equations)

                    @inbounds for v in eachvariable(equations)
                        du[v, i, j, element] += d_i_ii * fluxtilde1[v]
                        du[v, ii, j, element] += d_ii_i * fluxtilde1[v]
                    end
            end

            for jj in (j+1):nnodes(dg)
                d_j_jj = D[j, jj]
                d_jj_j = D[jj, j]

                if iszero(d_j_jj) && iszero(d_jj_j)
                    continue
                end

                u_node_jj = u_element[i, jj]
                Ja2_node_jj = get_contravariant_vector(2, contravariant_vectors, i, jj, element)
                Ja2_avg = (Ja2_node + Ja2_node_jj) * 0.5f0

                fluxtilde2 = volume_flux(u_node, u_node_jj, Ja2_avg, equations)

                    @inbounds for v in eachvariable(equations)
                        du[v, i, j, element] += d_j_jj * fluxtilde2[v]
                        du[v, i, jj, element] += d_jj_j * fluxtilde2[v]
                    end
            end
        end
    end

    return nothing
end

@inline function Trixi.calc_volume_integral!(du, u,
    mesh::UnstructuredMesh2D,
    have_nonconservative_terms::True, equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::FDSBP, cache)
    calc_volume_integral!(du, u, mesh, have_nonconservative_terms, combine_conservative_and_nonconservative_fluxes(volume_integral.volume_flux, equations),
        equations,
        volume_integral,
        dg, cache)

    return nothing
end

function calc_volume_integral!(du, u,
    mesh::UnstructuredMesh2D,
    have_nonconservative_terms::True, combine_conservative_and_nonconservative_fluxes::True, equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::FDSBP, cache)
    D = cache.derivative_split # SBP derivative operator
    @unpack contravariant_vectors = cache.elements
    @unpack volume_flux = volume_integral

    if nvariables(equations) == 1
        u_vectors = reshape(reinterpret(SVector{nvariables(equations),eltype(u)}, u),
            nnodes(dg), nnodes(dg), nelements(dg, cache))
    else
        u_vectors = reinterpret(reshape, SVector{nvariables(equations),eltype(u)}, u)
    end

    @threaded for element in eachelement(dg, cache)
        u_element = view(u_vectors, :, :, element)

        @inbounds for j in eachnode(dg), i in eachnode(dg)
            u_node = u_element[i, j]
            Ja1_node = get_contravariant_vector(1, contravariant_vectors, i, j, element)
            Ja2_node = get_contravariant_vector(2, contravariant_vectors, i, j, element)

            for ii in (i+1):nnodes(dg)
                d_i_ii = D[i, ii]
                d_ii_i = D[ii, i]

                if iszero(d_i_ii) && iszero(d_ii_i)
                    continue
                end

                u_node_ii = u_element[ii, j]
                Ja1_node_ii = get_contravariant_vector(1, contravariant_vectors, ii, j, element)
                Ja1_avg = (Ja1_node + Ja1_node_ii) * 0.5f0

                fluxtilde1_left, fluxtilde1_right = volume_flux(u_node, u_node_ii, Ja1_avg, equations)

                if !iszero(d_i_ii)
                    @inbounds for v in eachvariable(equations)
                        du[v, i, j, element] += d_i_ii * fluxtilde1_left[v]
                    end
                end

                if !iszero(d_ii_i)
                    @inbounds for v in eachvariable(equations)
                        du[v, ii, j, element] += d_ii_i * fluxtilde1_right[v]
                    end
                end
            end

            for jj in (j+1):nnodes(dg)
                d_j_jj = D[j, jj]
                d_jj_j = D[jj, j]

                if iszero(d_j_jj) && iszero(d_jj_j)
                    continue
                end

                u_node_jj = u_element[i, jj]
                Ja2_node_jj = get_contravariant_vector(2, contravariant_vectors, i, jj, element)
                Ja2_avg = (Ja2_node + Ja2_node_jj) * 0.5f0

                fluxtilde2_left, fluxtilde2_right = volume_flux(u_node, u_node_jj, Ja2_avg, equations)

                if !iszero(d_j_jj)
                    @inbounds for v in eachvariable(equations)
                        du[v, i, j, element] += d_j_jj * fluxtilde2_left[v]
                    end
                end

                if !iszero(d_jj_j)
                    @inbounds for v in eachvariable(equations)
                        du[v, i, jj, element] += d_jj_j * fluxtilde2_right[v]
                    end
                end
            end
        end
    end

    return nothing
end

function Trixi.create_cache(mesh::Union{TreeMesh{2},UnstructuredMesh2D}, equations,
    volume_integral::VolumeIntegralFluxDifferencing,
    dg::FDSBP, cache_containers, uEltype)
    prototype = Trixi.Array{SVector{nvariables(equations),uEltype},ndims(mesh)}(undef, ntuple(_ -> Trixi.nnodes(dg), ndims(mesh))...)
    f_threaded = [similar(prototype) for _ in 1:Threads.maxthreadid()]

    D = Matrix(dg.basis)
    M = mass_matrix(dg.basis)
    derivative_split = 2 .* D
    derivative_split[1, 1] += 1 / M[1, 1] # B[1, 1] = -1
    derivative_split[end, end] -= 1 / M[end, end] # B[end, end] = 1
    return (; f_threaded, derivative_split)
end


function Trixi.calc_surface_integral!(du, u, mesh::UnstructuredMesh2D,
    equations, surface_integral::SurfaceIntegralStrongForm,
    dg::DG, cache)
    inv_weight_left = inv(left_boundary_weight(dg.basis))
    inv_weight_right = inv(right_boundary_weight(dg.basis))
    @unpack normal_directions, surface_flux_values = cache.elements
    f_node = zero(eltype(u))
    @threaded for element in eachelement(dg, cache)
        for l in eachnode(dg)
            # surface at -x
            u_node = get_node_vars(u, equations, dg, 1, l, element)
            # compute internal flux in normal direction on side 4
            outward_direction = get_surface_normal(normal_directions, l, 4, element)
            f_num = get_node_vars(surface_flux_values, equations, dg, l, 4, element)
            multiply_add_to_node_vars!(du, inv_weight_left, f_num,
                equations, dg, 1, l, element)

            # surface at +x
            u_node = get_node_vars(u, equations, dg, nnodes(dg), l, element)
            # compute internal flux in normal direction on side 2
            outward_direction = get_surface_normal(normal_directions, l, 2, element)
            f_num = get_node_vars(surface_flux_values, equations, dg, l, 2, element)
            multiply_add_to_node_vars!(du, inv_weight_right, f_num,
                equations, dg, nnodes(dg), l, element)

            # surface at -y
            u_node = get_node_vars(u, equations, dg, l, 1, element)
            # compute internal flux in normal direction on side 1
            outward_direction = get_surface_normal(normal_directions, l, 1, element)
            f_num = get_node_vars(surface_flux_values, equations, dg, l, 1, element)
            multiply_add_to_node_vars!(du, inv_weight_left, f_num,
                equations, dg, l, 1, element)

            # surface at +y
            u_node = get_node_vars(u, equations, dg, l, nnodes(dg), element)
            # compute internal flux in normal direction on side 3
            outward_direction = get_surface_normal(normal_directions, l, 3, element)
            f_num = get_node_vars(surface_flux_values, equations, dg, l, 3, element)
            multiply_add_to_node_vars!(du, inv_weight_right, f_num,
                equations, dg, l, nnodes(dg), element)
        end
    end

    return nothing
end

function Trixi.calc_interface_flux!(surface_flux_values,
    mesh::UnstructuredMesh2D,
    have_nonconservative_terms::True, equations,
    surface_integral, dg::DG, cache)

    calc_interface_flux!(surface_flux_values, mesh, have_nonconservative_terms, combine_conservative_and_nonconservative_fluxes(surface_integral.surface_flux, equations), equations, surface_integral, dg, cache)

    return nothing
end



# compute the numerical flux interface with nonconservative terms coupling between two elements
# on an unstructured quadrilateral mesh
function calc_interface_flux!(surface_flux_values,
    mesh::UnstructuredMesh2D,
    have_nonconservative_terms::True, combine_conservative_and_nonconservative_fluxes::True,
    equations, surface_integral, dg::DG, cache)
    surface_flux = surface_integral.surface_flux
    @unpack u, start_index, index_increment, element_ids, element_side_ids = cache.interfaces
    @unpack normal_directions = cache.elements

    @threaded for interface in eachinterface(dg, cache)
        # Get the primary element index and local side index
        primary_element = element_ids[1, interface]
        primary_side = element_side_ids[1, interface]

        # Get neighboring element, local side index, and index increment on the
        # secondary element
        secondary_element = element_ids[2, interface]
        secondary_side = element_side_ids[2, interface]
        secondary_index_increment = index_increment[interface]

        secondary_index = start_index[interface]
        for primary_index in eachnode(dg)
            # pull the primary and secondary states from the boundary u values
            u_ll = get_one_sided_surface_node_vars(u, equations, dg, 1, primary_index,
                interface)
            u_rr = get_one_sided_surface_node_vars(u, equations, dg, 2, secondary_index,
                interface)

            # pull the outward pointing (normal) directional vector
            # Note! This assumes a conforming approximation, more must be done in terms
            # of the normals for hanging nodes and other non-conforming approximation spaces
            outward_direction = get_surface_normal(normal_directions, primary_index,
                primary_side, primary_element)

            # Calculate the conservative portion of the numerical flux
            # Call pointwise numerical flux with rotation. Direction is normalized
            # inside this function
            flux_left, flux_right = surface_flux(u_ll, u_rr, outward_direction, equations)

            # Copy flux to primary and secondary element storage
            # Note the sign change for the components in the secondary element!
            for v in eachvariable(equations)
                # Note the factor 0.5 necessary for the nonconservative fluxes based on
                # the interpretation of global SBP operators coupled discontinuously via
                # central fluxes/SAT
                surface_flux_values[v, primary_index, primary_side, primary_element] = flux_left[v]
                surface_flux_values[v, secondary_index, secondary_side, secondary_element] = -flux_right[v]
            end

            # increment the index of the coordinate system in the secondary element
            secondary_index += secondary_index_increment
        end
    end

    return nothing
end


@inline function Trixi.calc_boundary_flux!(surface_flux_values, t, boundary_condition,
                                     mesh::UnstructuredMesh2D,
                                     have_nonconservative_terms::True, equations,
                                     surface_integral, dg::DG, cache,
                                     node_index, side_index, element_index,
                                     boundary_index)
    @unpack normal_directions = cache.elements
    @unpack u, node_coordinates = cache.boundaries

    # pull the inner solution state from the boundary u values on the boundary element
    u_inner = get_node_vars(u, equations, dg, node_index, boundary_index)

    # pull the outward pointing (normal) directional vector
    outward_direction = get_surface_normal(normal_directions, node_index, side_index,
                                           element_index)

    # get the external solution values from the prescribed external state
    x = get_node_coords(node_coordinates, equations, dg, node_index, boundary_index)

    # Call pointwise numerical flux functions for the conservative and nonconservative part
    # in the normal direction on the boundary
    flux  = boundary_condition(u_inner, outward_direction, x, t,
                                            surface_integral.surface_flux, equations)

    for v in eachvariable(equations)
        # Note the factor 0.5 necessary for the nonconservative fluxes based on
        # the interpretation of global SBP operators coupled discontinuously via
        # central fluxes/SATs
        surface_flux_values[v, node_index, side_index, element_index] = flux[v]
    end

    return nothing
end
