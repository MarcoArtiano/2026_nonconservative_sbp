@inline function flux_zero(u_ll, u_rr, orientation, equations)
    zeroT = zero(eltype(u_ll))
    return zero(u_ll)
end