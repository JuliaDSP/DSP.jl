# This file was formerly a part of Julia. License is MIT: https://julialang.org/license

import Base: trailingsize

const SMALL_FILT_CUTOFF = 58

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
        elseif bs <= SMALL_FILT_CUTOFF
            _small_filt_fir!(out, b, x, si, col)
        else
            _filt_fir!(out, b, x, si, col)
        end
    end
    return out
end

# Transposed direct form II
function _filt_iir!(out, b, a, x, si, col)
    silen = length(si)
    @inbounds for i in axes(x, 1)
        xi = x[i,col]
        val = muladd(xi, b[1], si[1])
        for j=1:(silen-1)
            si[j] = muladd(val, -a[j+1], muladd(xi, b[j+1], si[j+1]))
        end
        si[silen] = muladd(xi, b[silen+1], -a[silen+1]*val)
        out[i,col] = val
    end
end

# Transposed direct form II
function _filt_fir!(out, b, x, si, col)
    silen = length(si)
    @inbounds for i in axes(x, 1)
        xi = x[i,col]
        val = muladd(xi, b[1], si[1])
        for j=1:(silen-1)
            si[j] = muladd(xi, b[j+1], si[j+1])
        end
        si[silen] = b[silen+1]*xi
        out[i,col] = val
    end
end

#
# filt implementation for FIR filters (faster than Base)
#

for n = 2:SMALL_FILT_CUTOFF
    silen = n-1
    si = [Symbol("si$i") for i = 1:silen]
    # Transposed direct form II
    @eval function _filt_fir!(out, b::NTuple{$n,T}, x, siarr, col) where {T}
        offset = (col - 1) * size(x, 1)

        $(Expr(:block, [:(@inbounds $(si[i]) = siarr[$i]) for i = 1:silen]...))
        @inbounds for i in axes(x, 1)
            xi = x[i+offset]
            val = muladd(xi, b[1], $(si[1]))
            $(Expr(:block, [:($(si[j]) = muladd(xi, b[$(j+1)], $(si[j+1]))) for j = 1:(silen-1)]...))
            $(si[silen]) = b[$(silen+1)]*xi
            out[i+offset] = val
        end
    end
end

# Convert array filter tap input to tuple for small-filtering
let chain = :(throw(ArgumentError("invalid tuple size")))
    for n = SMALL_FILT_CUTOFF:-1:2
        chain = quote
            if length(h) == $n
                _filt_fir!(
                    out,
                    ($([:(@inbounds(h[$i])) for i = 1:n]...),),
                    x,
                    si,
                    col
                )
            else
                $chain
            end
        end
    end

    @eval function _small_filt_fir!(
        out::AbstractArray, h::AbstractVector{T}, x::AbstractArray, si, col
    ) where T
        $chain
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


"""
    _zeropad!(padded::AbstractVector,
              u::AbstractVector,
              padded_axes = axes(padded),
              data_dest::Tuple = (1,),
              data_region = CartesianIndices(u))

Place the portion of `u` specified by `data_region` into `padded`, starting at
location `data_dest`. Sets the rest of `padded` to zero. This will mutate
`padded`. `padded_axes` must correspond to the axes of `padded`.

"""
@inline function _zeropad!(
    padded::AbstractVector,
    u::AbstractVector,
    padded_axes = axes(padded),
    data_dest::Tuple = (first(padded_axes[1]),),
    data_region = CartesianIndices(u),
)
    datasize = length(data_region)
    # Use axes to accommodate arrays that do not start at index 1
    data_first_i = first(data_region)[1]
    dest_first_i = data_dest[1]
    copyto!(padded, dest_first_i, u, data_first_i, datasize)
    padded[first(padded_axes[1]):dest_first_i - 1] .= 0
    padded[dest_first_i + datasize : end] .= 0

    padded
end

@inline function _zeropad!(
    padded::AbstractArray,
    u::AbstractArray,
    padded_axes = axes(padded),
    data_dest::Tuple = first.(padded_axes),
    data_region = CartesianIndices(u),
)
    fill!(padded, zero(eltype(padded)))
    dest_axes = UnitRange.(data_dest, data_dest .+ size(data_region) .- 1)
    dest_region = CartesianIndices(dest_axes)
    copyto!(padded, dest_region, u, data_region)

    padded
end

"""
    _zeropad(u, padded_size, [data_dest, data_region])

Creates and returns a new base-1 index array of size `padded_size`, with the
section of `u` specified by `data_region` copied into the region of the new
 array as specified by `data_dest`. All other values will be initialized to
 zero.

If either `data_dest` or `data_region` is not specified, then the defaults
described in [`_zeropad!`](@ref) will be used.
"""
function _zeropad(u, padded_size, args...)
    padded = similar(u, padded_size)
    _zeropad!(padded, u, axes(padded), args...)
end

function _zeropad_keep_offset(u, padded_size, u_axes, args...)
    ax_starts = first.(u_axes)
    new_axes = UnitRange.(ax_starts, ax_starts .+ padded_size .- 1)
    _zeropad!(similar(u, new_axes), u, args...)
end

function _zeropad_keep_offset(
    u, padded_size, ::NTuple{<:Any, Base.OneTo{Int}}, args...
)
    _zeropad(u, padded_size, args...)
end

"""
    _zeropad_keep_offset(u, padded_size, [data_dest, dat_region])

Like [`_zeropad`](@ref), but retains the first index of `u` when creating a new
array.
"""
function _zeropad_keep_offset(u, padded_size, args...)
    _zeropad_keep_offset(u, padded_size, axes(u), args...)
end

"""
Estimate the number of floating point multiplications per output sample for an
overlap-save algorithm with fft size `nfft`, and filter size `nb`.
"""
os_fft_complexity(nfft, nb) =  (nfft * log2(nfft) + nfft) / (nfft - nb + 1)

"""
Determine the length of FFT that minimizes the number of multiplications per
output sample for an overlap-save convolution of vectors of size `nb` and `nx`.
"""
function optimalfftfiltlength(nb, nx)
    nfull = nb + nx - 1

    # Step through possible nffts and find the nfft that minimizes complexity
    # Assumes that complexity is convex
    first_pow2 = ceil(Int, log2(nb))
    max_pow2 = ceil(Int, log2(nfull))
    prev_complexity = os_fft_complexity(2 ^ first_pow2, nb)
    pow2 = first_pow2 + 1
    while pow2 <= max_pow2
        new_complexity = os_fft_complexity(2 ^ pow2, nb)
        new_complexity > prev_complexity && break
        prev_complexity = new_complexity
        pow2 += 1
    end
    nfft = pow2 > max_pow2 ? 2 ^ max_pow2 : 2 ^ (pow2 - 1)

    if nfft > nfull
        # If nfft > nfull, then it's better to find next fast power
        nfft = nextfastfft(nfull)
    end

    nfft
end

"""
Prepare buffers and FFTW plans for convolution. The two buffers, tdbuff and
fdbuff may be an alias of each other, because complex convolution only requires
one buffer. The plans are mutating where possible, and the inverse plan is
unnormalized.
"""
@inline function os_prepare_conv(u::AbstractArray{T, N},
                                 nffts) where {T<:Real, N}
    tdbuff = similar(u, nffts)
    bufsize = ntuple(i -> i == 1 ? nffts[i] >> 1 + 1 : nffts[i], N)
    fdbuff = similar(u, Complex{T}, NTuple{N, Int}(bufsize))

    p = plan_rfft(tdbuff)
    ip = plan_brfft(fdbuff, nffts[1])

    tdbuff, fdbuff, p, ip
end

@inline function os_prepare_conv(u::AbstractArray{<:Complex}, nffts)
    buff = similar(u, nffts)

    p = plan_fft!(buff)
    ip = plan_bfft!(buff)

    buff, buff, p, ip # Only one buffer for complex
end

"""
Transform the smaller convolution input to frequency domain, and return it in a
new array. However, the contents of `buff` may be modified.
"""
@inline function os_filter_transform!(buff::AbstractArray{<:Real}, p)
    p * buff
end

@inline function os_filter_transform!(buff::AbstractArray{<:Complex}, p!)
    copy(p! * buff) # p operates in place on buff
end

"""
Take a block of data, and convolve it with the smaller convolution input. This
may modify the contents of `tdbuff` and `fdbuff`, and the result will be in
`tdbuff`.
"""
@inline function os_conv_block!(tdbuff::AbstractArray{<:Real},
                                fdbuff::AbstractArray,
                                filter_fd,
                                p,
                                ip)
    mul!(fdbuff, p, tdbuff)
    fdbuff .*= filter_fd
    mul!(tdbuff, ip, fdbuff)
end

"Like the real version, but only operates on one buffer"
@inline function os_conv_block!(buff::AbstractArray{<:Complex},
                                ::AbstractArray, # Only one buffer for complex
                                filter_fd,
                                p!,
                                ip!)
    p! * buff # p! operates in place on buff
    buff .*= filter_fd
    ip! * buff # ip! operates in place on buff
end

# Used by `unsafe_conv_kern_os!` to handle blocks of input data that need to be padded.
#
# For a given number of edge dimensions, convolve all regions along the
# perimeter that have that number of edge dimensions
#
# For a 3d cube, if n_edges = 1, this means all faces. If n_edges = 2, then
# all edges. Finally, if n_edges = 3, then all corners.
#
# This needs to be a separate function for subsets to generate tuple elements,
# which is only the case if the number of edges is a Val{n} type. Iterating over
# the number of edges with Val{n} is inherently type unstable, so this function
# boundary allows dispatch to make efficient code for each number of edge
# dimensions.
function unsafe_conv_kern_os_edge!(
    # These arrays and buffers will be mutated
    out::AbstractArray{<:Any, N},
    tdbuff,
    fdbuff,
    perimeter_range,
    # Number of edge dimensions to pad and convolve
    n_edges::Val,
    # Data to be convolved
    u,
    filter_fd,
    # FFTW plans
    p,
    ip,
    # Sizes, ranges, and other pre-calculated constants
    #
    ## ranges listing center and edge blocks
    edge_ranges,
    center_block_ranges,
    ## size and axis information
    all_dims, # 1:N
    su,
    u_start,
    sv,
    nffts,
    out_start,
    out_stop,
    save_blocksize,
    sout_deficit, # How many output samples are missing if nffts > sout
    tdbuff_axes,
) where N
    # Iterate over all possible combinations of edge dimensions for a number of
    # edges
    #
    # For a 3d cube with n_edges = 1, this will specify the top and bottom faces
    # (edge_dims = (1,)), then the left and right faces (edge_dims = (2,)), then
    # the front and back faces (edge_dims = (3,))
    for edge_dims in subsets(all_dims, n_edges)
        # Specify a region on the perimeter by combining an edge block index for
        # each dimension on an edge, and the central blocks for dimensions not
        # on an edge.
        #
        # First make all entries equal to the center blocks:
        @inbounds copyto!(perimeter_range, 1, center_block_ranges, 1, N)

        # For the dimensions chosen to be on an edge (edge_dims), get the
        # ranges of the blocks that would need to be padded (lie on an edge)
        # in that dimension.
        #
        # There can be one to two such ranges for each dimension, because with
        # some inputs sizes the whole convolution is just one range
        # (one edge block), or the padding will be needed at both the leading
        # and trailing side of that dimension
        selected_edge_ranges = getindex.(Ref(edge_ranges), edge_dims)

        # Visit each combination of edge ranges for the edge dimensions chosen.
        # For a 3d cube with n_edges = 1 and edge_dims = (1,), this will visit
        # the top face, and then the bottom face.
        for perimeter_edge_ranges in Iterators.ProductIterator(selected_edge_ranges)
            # The center region for non-edge dimensions has been specified above,
            # so finish specifying the region of the perimeter for this edge
            # block
            @inbounds for (i, dim) in enumerate(edge_dims)
                perimeter_range[dim] = perimeter_edge_ranges[i]
            end

            # Region of block indices, not data indices!
            block_region = CartesianIndices(
                NTuple{N, UnitRange{Int}}(perimeter_range)
            )
            @inbounds for block_pos in block_region
                # Figure out which portion of the input data should be transformed

                block_idx = convert(NTuple{N, Int}, block_pos)
                ## data_offset is NOT the beginning of the region that will be
                ## convolved, but is instead the beginning of the unaliased data.
                data_offset = save_blocksize .* (block_idx .- 1)
                ## How many zeros will need to be added before the data
                pad_before = max.(0, sv .- data_offset .- 1)
                data_ideal_stop = data_offset .+ save_blocksize
                ## How many zeros will need to be added after the data
                pad_after = max.(0, data_ideal_stop .- su)

                ## Data indices, not block indices
                data_region = CartesianIndices(
                    UnitRange.(
                        u_start .+ data_offset .- sv .+ pad_before .+ 1,
                        u_start .+ data_ideal_stop .- pad_after .- 1
                    )
                )

                # Convolve portion of input

                _zeropad!(tdbuff, u, tdbuff_axes, pad_before .+ 1, data_region)
                os_conv_block!(tdbuff, fdbuff, filter_fd, p, ip)

                # Save convolved result to output

                block_out_stop = min.(
                    out_start .+ data_offset .+ save_blocksize .- 1,
                    out_stop
                )
                block_out_region = CartesianIndices(
                    UnitRange.(out_start .+ data_offset, block_out_stop)
                )
                ## If the input could not fill tdbuff, account for that before
                ## copying the convolution result to the output
                u_deficit = max.(0, pad_after .- sv .+ 1)
                valid_buff_region = CartesianIndices(
                    UnitRange.(sv, nffts .- u_deficit .- sout_deficit)
                )
                copyto!(out, block_out_region, tdbuff, valid_buff_region)
            end
        end
    end
end

# Assumes u is larger than, or the same size as, v
# nfft should be greater than or equal to 2*sv-1
function unsafe_conv_kern_os!(out,
                        u::AbstractArray{<:Any, N},
                        v,
                        su,
                        sv,
                        sout,
                        nffts) where N
    u_start = first.(axes(u))
    out_axes = axes(out)
    out_start = first.(out_axes)
    out_stop = last.(out_axes)
    ideal_save_blocksize = nffts .- sv .+ 1
    # Number of samples that are "missing" if the output is smaller than the
    # valid portion of the convolution
    sout_deficit = max.(0, ideal_save_blocksize .- sout)
    # Size of the valid portion of the convolution result
    save_blocksize = ideal_save_blocksize .- sout_deficit
    nblocks = cld.(sout, save_blocksize)

    # Pre-allocation
    tdbuff, fdbuff, p, ip = os_prepare_conv(u, nffts)
    tdbuff_axes = axes(tdbuff)

    # Transform the smaller filter
    _zeropad!(tdbuff, v)
    filter_fd = os_filter_transform!(tdbuff, p)
    filter_fd .*= 1 / prod(nffts) # Normalize once for brfft

    # block indices for center blocks, which need no padding
    first_center_blocks = cld.(sv .- 1, save_blocksize) .+ 1
    last_center_blocks = fld.(su, save_blocksize)
    center_block_ranges = UnitRange.(first_center_blocks, last_center_blocks)

    # block index ranges for blocks that need to be padded
    # Corresponds to the leading and trailing side of a dimension, or if there
    # are no center blocks, corresponds to the whole dimension
    edge_ranges = map(nblocks, first_center_blocks, last_center_blocks) do nblock, firstfull, lastfull
        lastfull > 1 ? [1:firstfull - 1, lastfull + 1 : nblock] : [1:nblock]
    end
    all_dims = 1:N
    # Buffer to store ranges of indices for a single region of the perimeter
    perimeter_range = Vector{UnitRange{Int}}(undef, N)

    # Convolve all blocks that require padding.
    #
    # This is accomplished by dividing the perimeter of the volume into
    # subsections, where the division is done by the number of edge dimensions.
    # For a 3d cube, this convolves the perimeter in the following order:
    #
    # Number of Edge Dimensions | Convolved Region
    # --------------------------+-----------------
    #                         1 | Faces of Cube
    #                         2 | Edges of Cube
    #                         3 | Corners of Cube
    #
    for n_edges in all_dims
        unsafe_conv_kern_os_edge!(
            # These arrays and buffers will be mutated
            out,
            tdbuff,
            fdbuff,
            perimeter_range,
            # Number of edge dimensions to pad and convolve
            Val(n_edges),
            # Data to be convolved
            u,
            filter_fd,
            # FFTW plans
            p,
            ip,
            # Sizes, ranges, and other pre-calculated constants
            #
            ## ranges listing center and edge blocks
            edge_ranges,
            center_block_ranges,
            ## size and axis information
            all_dims, # 1:N
            su,
            u_start,
            sv,
            nffts,
            out_start,
            out_stop,
            save_blocksize,
            sout_deficit,
            tdbuff_axes) # How many output samples are missing if nffts > sout
    end

    tdbuff_region = CartesianIndices(tdbuff)
    # Portion of buffer with valid result of convolution
    valid_buff_region = CartesianIndices(UnitRange.(sv, nffts))
    # Iterate over block indices (not data indices) that do not need to be padded
    @inbounds for block_pos in CartesianIndices(center_block_ranges)
        # Calculate portion of data to transform

        block_idx = convert(NTuple{N, Int}, block_pos)
        ## data_offset is NOT the beginning of the region that will be
        ## convolved, but is instead the beginning of the unaliased data.
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

function _conv_kern_fft!(out,
                         u::AbstractArray{T, N},
                         v::AbstractArray{T, N},
                         su,
                         sv,
                         outsize,
                         nffts) where {T<:Real, N}
    padded = _zeropad(u, nffts)
    p = plan_rfft(padded)
    uf = p * padded
    _zeropad!(padded, v)
    vf = p * padded
    uf .*= vf
    raw_out = irfft(uf, nffts[1])
    copyto!(out,
            CartesianIndices(out),
            raw_out,
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
    copyto!(out,
            CartesianIndices(out),
            upad,
            CartesianIndices(UnitRange.(1, outsize)))
end

# v should be smaller than u for good performance
function _conv_fft!(out, u, v, su, sv, outsize)
    os_nffts = map(optimalfftfiltlength, sv, su)
    if any(os_nffts .< outsize)
        unsafe_conv_kern_os!(out, u, v, su, sv, outsize, os_nffts)
    else
        nffts = nextfastfft(outsize)
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

# We use this type definition for clarity
const RealOrComplexFloat = Union{AbstractFloat, Complex{T} where T<:AbstractFloat}

# May switch argument order
"""
    conv(u,v)

Convolution of two arrays. Uses either FFT convolution or overlap-save,
depending on the size of the input. `u` and `v` can be  N-dimensional arrays,
with arbitrary indexing offsets, but their axes must be a `UnitRange`.
"""
function conv(u::AbstractArray{T, N},
              v::AbstractArray{T, N}) where {T<:RealOrComplexFloat, N}
    su = size(u)
    sv = size(v)
    if prod(su) >= prod(sv)
        _conv(u, v, su, sv)
    else
        _conv(v, u, sv, su)
    end
end

function conv(u::AbstractArray{<:RealOrComplexFloat, N},
              v::AbstractArray{<:RealOrComplexFloat, N}) where N
    fu, fv = promote(u, v)
    conv(fu, fv)
end

conv(u::AbstractArray{<:Integer, N}, v::AbstractArray{<:Integer, N}) where {N} =
    round.(Int, conv(float(u), float(v)))

conv(u::AbstractArray{<:Number, N}, v::AbstractArray{<:Number, N}) where {N} =
    conv(float(u), float(v))

function conv(u::AbstractArray{<:Number, N},
              v::AbstractArray{<:RealOrComplexFloat, N}) where N
    conv(float(u), v)
end

function conv(u::AbstractArray{<:RealOrComplexFloat, N},
              v::AbstractArray{<:Number, N}) where N
    conv(u, float(v))
end

function conv(A::AbstractArray{<:Number, M},
              B::AbstractArray{<:Number, N}) where {M, N}
    if (M < N)
        conv(cat(A, dims=N)::AbstractArray{eltype(A), N}, B)
    else
        @assert M > N
        conv(A, cat(B, dims=M)::AbstractArray{eltype(B), M})
    end
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


dsp_reverse(v, ::NTuple{<:Any, Base.OneTo{Int}}) = reverse(v, dims = 1)
function dsp_reverse(v, vaxes)
    vsize = length(v)
    reflected_start = - first(vaxes[1]) - vsize + 1
    reflected_axes = (reflected_start : reflected_start + vsize - 1,)
    out = similar(v, reflected_axes)
    copyto!(out, reflected_start, Iterators.reverse(v), 1, vsize)
end


"""
    xcorr(u; padmode::Symbol=:none, scaling::Symbol=:none)
    xcorr(u, v; padmode::Symbol=:none, scaling::Symbol=:none)

With two arguments, compute the cross-correlation of two vectors, by calculating
the similarity between `u` and `v` with various offsets of `v`. Delaying `u`
relative to `v` will shift the result to the right. If one argument is provided,
calculate `xcorr(u, u; kwargs...)`.

The size of the output depends on the `padmode` keyword argument: with `padmode =
:none` the length of the result will be `length(u) + length(v) - 1`, as with `conv`.
With `padmode = :longest`, the shorter of the arguments will be padded so they
are of equal length. This gives a result with length `2*max(length(u), length(v))-1`,
with the zero-lag condition at the center.

The keyword argument `scaling` can be provided. Possible arguments are the default
`:none` and `:biased`. `:biased` is valid only if the vectors have the same length,
or only one vector is provided, dividing the result by `length(u)`.
"""
function xcorr(
    u::AbstractVector, v::AbstractVector; padmode::Symbol=:none, scaling::Symbol=:none
)
    su = size(u, 1); sv = size(v, 1)
    res = if padmode == :longest
        if su < sv
            u = _zeropad_keep_offset(u, sv)
        elseif sv < su
            v = _zeropad_keep_offset(v, su)
        end
        conv(u, dsp_reverse(conj(v), axes(v)))
    elseif padmode == :none
        conv(u, dsp_reverse(conj(v), axes(v)))
    else
        throw(ArgumentError("padmode keyword argument must be either :none or :longest"))
    end

    if scaling == :biased
        if su == sv
            res ./= su
        else
            throw(DimensionMismatch("scaling only valid for vectors of same length"))
        end
    end

    res
end

xcorr(u::AbstractVector; kwargs...) = xcorr(u, u; kwargs...)
