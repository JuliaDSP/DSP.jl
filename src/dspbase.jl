# This file was formerly a part of Julia. License is MIT: https://julialang.org/license

import Base.trailingsize
import LinearAlgebra.BLAS

_zerosi(b,a,T) = zeros(promote_type(eltype(b), eltype(a), T), max(length(a), length(b))-1)

"""
    filt(b, a, x, [si])

Apply filter described by vectors `a` and `b` to vector `x`, with an optional initial filter
state vector `si` (defaults to zeros).
"""
function filt(b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number},
              x::AbstractArray{T}, si::AbstractArray{S} = _zerosi(b,a,T)) where {T,S}
    filt!(Array{promote_type(eltype(b), eltype(a), T, S)}(undef, size(x)), b, a, x, si)
end

# in-place filtering: returns results in the out argument, which may shadow x
# (and does so by default)

"""
    filt!(out, b, a, x, [si])

Same as [`filt`](@ref) but writes the result into the `out` argument, which may
alias the input `x` to modify it in-place.
"""
function filt!(out::AbstractArray, b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number},
               x::AbstractArray{T}, si::AbstractArray{S,N} = _zerosi(b,a,T)) where {T,S,N}
    isempty(b) && throw(ArgumentError("filter vector b must be non-empty"))
    isempty(a) && throw(ArgumentError("filter vector a must be non-empty"))
    a[1] == 0  && throw(ArgumentError("filter vector a[1] must be nonzero"))
    if size(x) != size(out)
        throw(ArgumentError("output size $(size(out)) must match input size $(size(x))"))
    end

    as = length(a)
    bs = length(b)
    sz = max(as, bs)
    silen = sz - 1
    ncols = trailingsize(x,2)

    if size(si, 1) != silen
        throw(ArgumentError("initial state vector si must have max(length(a),length(b))-1 rows"))
    end
    if N > 1 && trailingsize(si,2) != ncols
        throw(ArgumentError("initial state vector si must be a vector or have the same number of columns as x"))
    end

    size(x,1) == 0 && return out
    sz == 1 && return mul!(out, x, b[1]/a[1]) # Simple scaling without memory

    # Filter coefficient normalization
    if a[1] != 1
        norml = a[1]
        a = a ./ norml
        b = b ./ norml
    end
    # Pad the coefficients with zeros if needed
    bs<sz   && (b = copyto!(zeros(eltype(b), sz), b))
    1<as<sz && (a = copyto!(zeros(eltype(a), sz), a))

    initial_si = si
    for col = 1:ncols
        # Reset the filter state
        si = initial_si[:, N > 1 ? col : 1]
        if as > 1
            _filt_iir!(out, b, a, x, si, col)
        else
            _filt_fir!(out, b, x, si, col)
        end
    end
    return out
end

function _filt_iir!(out, b, a, x, si, col)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,col]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi - a[j+1]*val
        end
        si[silen] = b[silen+1]*xi - a[silen+1]*val
        out[i,col] = val
    end
end

function _filt_fir!(out, b, x, si, col)
    silen = length(si)
    @inbounds for i=1:size(x, 1)
        xi = x[i,col]
        val = si[1] + b[1]*xi
        for j=1:(silen-1)
            si[j] = si[j+1] + b[j+1]*xi
        end
        si[silen] = b[silen+1]*xi
        out[i,col] = val
    end
end

"""
    deconv(b,a) -> c

Construct vector `c` such that `b = conv(a,c) + r`.
Equivalent to polynomial division.
"""
function deconv(b::StridedVector{T}, a::StridedVector{T}) where T
    lb = size(b,1)
    la = size(a,1)
    if lb < la
        return [zero(T)]
    end
    lx = lb-la+1
    x = zeros(T, lx)
    x[1] = 1
    filt(b, a, x)
end

last_aligned_location(padded, u) = size(padded) .- size(u) .+ 1

# padded must start at index 1, but u can have arbitrary offset
function _zeropad!(
    padded::AbstractVector,
    u::AbstractVector,
    data_dest = (1,),
    data_region = CartesianIndices(u)
)
    datasize = length(data_region)
    start_i = data_dest[1]

    # Use axes to accommodate arrays that do not start at index 1
    data_first_i = first(data_region)[1]
    copyto!(padded, start_i, u, data_first_i, datasize)

    padded[1:start_i - 1] .= 0
    padded[start_i + datasize:end] .= 0

    padded
end

# padded must start at index 1, but u can have arbitrary offset
function _zeropad!(
    padded::AbstractArray{<:Any, N}, u::AbstractArray{<:Any, N},
    data_dest = ntuple(::Integer -> 1, N),
    data_region = CartesianIndices(u)
) where N
    # Copy the data to the beginning of the padded array
    fill!(padded, zero(eltype(padded)))
    pad_dest_axes = UnitRange.(data_dest, data_dest .+ size(data_region) .- 1)
    copyto!(padded, CartesianIndices(pad_dest_axes), u, data_region)

    padded
end
_zeropad(u, padded_size, args...) = _zeropad!(similar(u, padded_size), u, args...)

# Number of real operations required for overlap-save with nfft = 2^pow2 and filter
# length nb
os_fft_complexity(pow2, nb) = 4 * (2 ^ pow2 * (pow2 + 1)) / (2 ^ pow2 - nb + 1)

# Determine optimal length of the FFT for fftfilt
@inline function optimalfftfiltlength(nb, nx)
    nu, nv = ifelse(nb <= nx, (nb, nx), (nx, nb))
    first_pow2 = ceil(Int, log2(nu))
    last_pow2 = ceil(Int, log2(nv + nu - 1))
    complexities = os_fft_complexity.(first_pow2:last_pow2, nu)

    # Find power of 2 with least complexity relative to the first power of 2
    relative_ind_best_pow2 = argmin(complexities)

    best_pow2 = first_pow2 + relative_ind_best_pow2 - 1
    nfft = 2 ^ best_pow2

    L = nfft - nu + 1
    if L > nv
        # If L > nv, better to find next fast power
        nfft = nextfastfft(nv + nu - 1)
    end

    nfft
end

@inline function os_prepare_conv(u::AbstractArray{T, N}, nffts) where {T<:Real, N}
    tdbuff = similar(u, nffts)
    size_arr = collect(nffts)
    size_arr[1] = size_arr[1] >> 1 + 1
    fdbuff = similar(u, Complex{T}, NTuple{N, Int}(size_arr))

    p = plan_rfft(tdbuff)
    ip = plan_brfft(fdbuff, nffts[1])

    tdbuff, fdbuff, p, ip
end
@inline function os_prepare_conv(u::AbstractArray{<:Complex, <:Any}, nffts)
    buff = similar(u, nffts)

    p = plan_fft!(buff)
    ip = plan_bfft!(buff)

    buff, buff, p, ip
end

@inline function os_filter_transform!(buff::AbstractArray{<:Real, <:Any}, p)
    p * buff
end
@inline function os_filter_transform!(buff::AbstractArray{<:Complex, <:Any}, p!)
    copy(p! * buff) # p operates in place on buff
end

@inline function os_conv_block!(
    tdbuff::AbstractArray{<:Real, <:Any}, fdbuff::AbstractArray, filter_fd, p, ip
)
    mul!(fdbuff, p, tdbuff)
    fdbuff .*= filter_fd
    mul!(tdbuff, ip, fdbuff)
end
@inline function os_conv_block!(
    buff::AbstractArray{<:Complex, <:Any}, ::AbstractArray, filter_fd, p!, ip!
)
    p! * buff
    buff .*= filter_fd
    ip! * buff
end

function _conv_kern_os_edge!(
    out::AbstractArray{<:Any, N},
    tdbuff,
    fdbuff,
    edge_range,
    p,
    ip,
    n_edges,
    u,
    filter_fd,
    center_block_ranges,
    edge_blocks,
    all_dims,
    save_blocksize,
    su,
    sv,
    nffts,
    out_start,
    out_stop,
    u_start
) where N
    for edge_dims in subsets(all_dims, n_edges)
        center_dims = setdiff(all_dims, edge_dims)
        edge_dims_arr = collect(edge_dims)
        for dim in center_dims
            edge_range[dim] = center_block_ranges[dim]
        end
        these_blocks = getindex.(Ref(edge_blocks), edge_dims)
        for edge_block_nos in Iterators.ProductIterator(these_blocks)
            edge_range[edge_dims_arr] .= UnitRange.(
                edge_block_nos, edge_block_nos
            )
            block_region = CartesianIndices(
                NTuple{N, UnitRange{Int}}(edge_range)
            )
            for block_pos in block_region
                # Figure out which portion of the input needs to be transformed
                block_idx = convert(NTuple{N, Int}, block_pos)
                data_offset = save_blocksize .* (block_idx .- 1)
                pad_before = max.(0, sv .- data_offset .- 1)
                data_ideal_stop = data_offset .+ save_blocksize
                pad_after = max.(0, data_ideal_stop .- su)
                data_region = CartesianIndices(
                    UnitRange.(
                        u_start .+ data_offset .- sv .+ pad_before .+ 1,
                        u_start .+ data_ideal_stop .- pad_after .- 1
                    )
                )

                # Convolve portion of input
                _zeropad!(tdbuff, u, pad_before .+ 1, data_region)
                os_conv_block!(tdbuff, fdbuff, filter_fd, p, ip)

                # Save convolved result to output
                block_out_stop = min.(
                    out_start .+ data_offset .+ save_blocksize .- 1,
                    out_stop
                )
                block_out_region = CartesianIndices(
                    UnitRange.(out_start .+ data_offset, block_out_stop)
                )
                deficit = max.(0, pad_after .- sv .+ 1)
                valid_buff_region = CartesianIndices(
                    UnitRange.(sv, nffts .- deficit)
                )
                copyto!(out, block_out_region, tdbuff, valid_buff_region)
            end
        end
    end
end

# Assumes u is larger, or same size, as v
function _conv_kern_os!(out, u::AbstractArray{<:Any, N}, v, su, sv, sout, nffts) where N
    u_start = first.(axes(u))
    out_axes = axes(out)
    out_start = first.(out_axes)
    out_stop = last.(out_axes)
    save_blocksize = min.(sout, nffts .- sv .+ 1)
    nblocks = cld.(sout, save_blocksize)

    # Pre-allocation
    tdbuff, fdbuff, p, ip = os_prepare_conv(u, nffts)

    # Transform the smaller filter
    _zeropad!(tdbuff, v)
    filter_fd = os_filter_transform!(tdbuff, p)
    filter_fd .*= 1 / prod(nffts) # Normalize once for brfft

    center_block_ranges = UnitRange.(2, nblocks .- 1)
    edge_blocks = map(nblock -> nblock > 1 ? [1, nblock] : [1], nblocks)
    all_dims = 1:N
    edge_range = Vector{UnitRange}(undef, N)
    for n_edges in all_dims
        _conv_kern_os_edge!(
            out,
            tdbuff,
            fdbuff,
            edge_range,
            p,
            ip,
            Val{n_edges}(),
            u,
            filter_fd,
            center_block_ranges,
            edge_blocks,
            all_dims,
            save_blocksize,
            su,
            sv,
            nffts,
            out_start,
            out_stop,
            u_start
        )
    end

    tdbuff_region = CartesianIndices(tdbuff)
    valid_buff_region = CartesianIndices(UnitRange.(sv, nffts))
    for block_pos in CartesianIndices(center_block_ranges)
        # Calculate portion of data to transform
        block_idx = convert(NTuple{N, Int}, block_pos)
        data_offset = save_blocksize .* (block_idx .- 1)
        data_stop = data_offset .+ save_blocksize
        data_region = CartesianIndices(
            UnitRange.(u_start .+ data_offset .- sv .+ 1, u_start .+ data_stop .- 1)
        )

        # Convolve this portion of the data
        copyto!(tdbuff, tdbuff_region, u, data_region)
        os_conv_block!(tdbuff, fdbuff, filter_fd, p, ip)

        # Save convolved result to output
        block_out_region = CartesianIndices(
            UnitRange.(data_offset .+ out_start, data_stop .+ out_start .- 1)
        )
        copyto!(out, block_out_region, tdbuff, valid_buff_region)
    end

    out
end

function _conv_kern_fft!(
    out, u::AbstractArray{T, N}, v::AbstractArray{T, N}, su, sv, outsize, nffts
) where {T<:Real, N}
    padded = _zeropad(u, nffts)
    p = plan_rfft(padded)
    uf = p * padded
    _zeropad!(padded, v)
    vf = p * padded
    uf .*= vf
    raw_out = irfft(uf, nffts[1])
    copyto!(out, CartesianIndices(out), raw_out,
            CartesianIndices(UnitRange.(1, outsize)))
end
function _conv_kern_fft!(out, u, v, su, sv, outsize, nffts)
    upad = _zeropad(u, nffts)
    vpad = _zeropad(v, nffts)
    p! = plan_fft!(upad)
    p! * upad # Operates in place on upad
    p! * vpad
    upad .*= vpad
    ifft!(upad)
    copyto!(out, CartesianIndices(out), upad,
            CartesianIndices(UnitRange.(1, outsize)))
end

function _conv_fft!(out, u, v, su, sv, outsize)
    nffts = map((nu, nv)-> optimalfftfiltlength(nu, nv), su, sv)
    if any(nffts .< su .+ sv .- 1) 
        _conv_kern_os!(out, u, v, su, sv, outsize, nffts)
    else
        _conv_kern_fft!(out, u, v, su, sv, outsize, nffts)
    end
end


# For arrays with weird offsets
function _conv_similar(u, outsize, axesu, axesv)
    out_offsets = first.(axesu) .+ first.(axesv)
    out_axes = UnitRange.(out_offsets, out_offsets .+ outsize .- 1)
    similar(u, out_axes)
end
function _conv_similar(
    u, outsize, ::NTuple{<:Any, Base.OneTo{Int}}, ::NTuple{<:Any, Base.OneTo{Int}}
)
    similar(u, outsize)
end
_conv_similar(u, v, outsize) = _conv_similar(u, outsize, axes(u), axes(v))

# Does convolution, will not switch argument order
function _conv!(out, u, v, su, sv, outsize)
    # TODO: Add spatial / time domain algorithm
    _conv_fft!(out, u, v, su, sv, outsize)
end

# Does convolution, will not switch argument order
function _conv(u, v, su, sv)
    outsize = su .+ sv .- 1
    out = _conv_similar(u, v, outsize)
    _conv!(out, u, v, su, sv, outsize)
end

# May switch argument order
"""
    conv(u,v)

Convolution of two arrays. Uses FFT algorithm.
"""
function conv(u::AbstractArray{T, N},
              v::AbstractArray{T, N}) where {T<:BLAS.BlasFloat, N}
    su = size(u)
    sv = size(v)
    if prod(su) >= prod(sv)
        _conv(u, v, su, sv)
    else
        _conv(v, u, sv, su)
    end
end

function conv(u::AbstractArray{<:BLAS.BlasFloat, N},
              v::AbstractArray{<:BLAS.BlasFloat, N}) where N
    fu, fv = promote(u, v)
    conv(fu, fv)
end
conv(u::AbstractArray{T, N}, v::AbstractArray{T, N}) where {T<:Number, N} =
    conv(float(u), float(v))
conv(u::AbstractArray{<:Integer, N}, v::AbstractArray{<:Integer, N}) where {N} =
    round.(Int, conv(float(u), float(v)))
function conv(u::AbstractArray{<:Number, N},
              v::AbstractArray{<:BLAS.BlasFloat, N}) where N
    conv(float(u), v)
end
function conv(u::AbstractArray{<:BLAS.BlasFloat, N},
              v::AbstractArray{<:Number, N}) where N
    conv(u, float(v))
end

function conv(A::AbstractArray, B::AbstractArray)
    maxnd = max(ndims(A), ndims(B))
    return conv(cat(A, dims=maxnd), cat(B, dims=maxnd))
end

"""
    conv(u,v,A)

2-D convolution of the matrix `A` with the 2-D separable kernel generated by
the vectors `u` and `v`.
Uses 2-D FFT algorithm.
"""
function conv(u::AbstractVector{T}, v::AbstractVector{T}, A::AbstractMatrix{T}) where T
    # Arbitrary indexing offsets not implemented
    @assert !Base.has_offset_axes(u, v, A)
    m = length(u)+size(A,1)-1
    n = length(v)+size(A,2)-1
    B = zeros(T, m, n)
    B[1:size(A,1),1:size(A,2)] = A
    u = fft([u;zeros(T,m-length(u))])
    v = fft([v;zeros(T,n-length(v))])
    C = ifft(fft(B) .* (u * transpose(v)))
    if T <: Real
        return real(C)
    end
    return C
end


function check_padmode_kwarg(padmode::Symbol, su::Integer, sv::Integer)
    if padmode == :default_longest
        if su != sv
            Base.depwarn(
            """
            The default value of `padmode` will be changing from `:longest` to
            `:none` in a future release of DSP. In preparation for this change,
            leaving `padmode` unspecified is currently deprecated. To keep
            current behavior specify `padmode=:longest`. To avoid this warning,
            specify padmode = :none or padmode = :longest where appropriate.
            """
                ,
                :xcorr
            )
        end
        :longest
    else
        padmode
    end
end

dsp_reverse(v, ::NTuple{<:Any, Base.OneTo{Int}}) = reverse(v, dims = 1)
function dsp_reverse(v, vaxes)
    vsize = length(v)
    reflected_start = - first(vaxes[1]) - vsize + 1
    reflected_axes = (reflected_start : reflected_start + vsize - 1,)
    out = similar(v, reflected_axes)
    copyto!(out, reflected_start, Iterators.reverse(v), 1, vsize)
end


"""
    xcorr(u,v; padmode = :longest)

Compute the cross-correlation of two vectors, by calculating the similarity
between `u` and `v` with various offsets of `v`. Delaying `u` relative to `v`
will shift the result to the right.

The size of the output depends on the padmode keyword argument: with padmode =
:none the length of the result will be length(u) + length(v) - 1, as with conv.
With padmode = :longest the shorter of the arguments will be padded so they are
equal length. This gives a result with length 2*max(length(u), length(v))-1,
with the zero-lag condition at the center.

!!! warning
    The default value of `padmode` will be changing from `:longest` to `:none`
    in a future release of DSP. In preparation for this change, leaving
    `padmode` unspecified is currently deprecated.
"""
function xcorr(
    u::AbstractVector, v::AbstractVector; padmode::Symbol = :default_longest
)
    su = size(u,1); sv = size(v,1)
    padmode = check_padmode_kwarg(padmode, su, sv)
    if padmode == :longest
        if su < sv
            u = _zeropad(u, sv)
        elseif sv < su
            v = _zeropad(v, su)
        end
        conv(u, dsp_reverse(conj(v), axes(v)))
    elseif padmode == :none
        conv(u, dsp_reverse(conj(v), axes(v)))
    else
        throw(ArgumentError("padmode keyword argument must be either :none or :longest"))
    end
end
