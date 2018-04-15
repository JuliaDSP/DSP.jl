module Unwrap

export unwrap, unwrap!

"""
    unwrap!(m; range=2pi, wrap_around=tuple(fill(false, N)...), seed=-1)

In-place version of unwrap(m; range, wrap_around, seed)
"""
unwrap!

"""
    unwrap(m::AbstractArray{T, N}; range=2pi, wrap_around=tuple(fill(false, N)...), seed=-1)

Assumes `m` to be an array of values of arbitrary dimension that
have been wrapped to be inside the given range (centered around
zero), and undoes the wrapping by identifying discontinuities. To operate on a
subset of dimensions of an array, simply pass in the corresponding `view`.

A common usage is for a phase measurement over time, such as when
comparing successive frames of a short-time-fourier-transform, as
each frame is wrapped to stay within (-pi, pi].

# Arguments
- `A::AbstractArray{T, N}`: Array to unwrap
- `range::Number`: Range of wrapped array
- `wrap_around::NTuple{N, Bool}`:  When an element of this tuple is `true`, the
    unwrapping process will consider the edges along the corresponding axis
    of the image to be connected
- `seed::Int`: Unwrapping of arrays with dimension > 1 uses a random initialization. This
    sets the seed of the RNG.
"""
unwrap

function unwrap(m::AbstractArray{T, N}, args...; kwargs...) where {T<:AbstractFloat, N}
    unwrap!(copy(m), args...; kwargs...)
end

function unwrap!(A::AbstractVector{T};
                 range::Number=2*convert(T, pi),
                 wrap_around::NTuple{1, Bool}=(false,),
                 seed::Int=-1) where T
    @inbounds previous_element = A[1]
    difference = zero(previous_element)
    periods = 0
    @inbounds for i in 2:length(A)
        difference = A[i] - previous_element
        periods += ifelse(difference >  range/2, -1, 0)
        periods += ifelse(difference < -range/2,  1, 0)
        previous_element = A[i]
        A[i] = previous_element + one(previous_element) * periods * range
    end
    return A
end

#= Algorithm based off of
 M. A. Herráez, D. R. Burton, M. J. Lalor, and M. A. Gdeisat,
 "Fast two-dimensional phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path"
 `Applied Optics, Vol. 41, Issue 35, pp. 7437-7444 (2002) <http://dx.doi.org/10.1364/AO.41.007437>`
 and
 H. Abdul-Rahman, M. Gdeisat, D. Burton, M. Lalor,
 "Fast three-dimensional phase-unwrapping algorithm based on sorting by reliability following a non-continuous path",
 `Proc. SPIE 5856, Optical Measurement Systems for Industrial Inspection IV, 32 (2005) <http://dx.doi.ogr/doi:10.1117/12.611415>`
 Code inspired by Scipy's implementation, which is under BSD license.
=#

mutable struct Pixel{T}
    periods::Int
    val::T
    reliability::Float64
    groupsize::Int
    head::Pixel{T}
    last::Pixel{T}
    next::Union{Void, Pixel{T}}
    function Pixel{T}(periods, val, rel, gs) where T
        pixel = new(periods, val, rel, gs)
        pixel.head = pixel
        pixel.last = pixel
        pixel.next = nothing
        return pixel
    end
end
Pixel(v) = Pixel{typeof(v)}(0, v, rand(), 1)
@inline Base.length(p::Pixel) = p.head.groupsize

struct Edge{N}
    reliability::Float64
    periods::Int
    pixel_1::CartesianIndex{N}
    pixel_2::CartesianIndex{N}
end
function Edge{N}(pixel_image::AbstractArray, ind1::CartesianIndex{N}, ind2::CartesianIndex{N}, range) where N
    @inbounds rel = pixel_image[ind1].reliability + pixel_image[ind2].reliability
    @inbounds periods = find_period(pixel_image[ind1].val, pixel_image[ind2].val, range)
    return Edge{N}(rel, periods, ind1, ind2)
end
@inline Base.isless(e1::Edge, e2::Edge) = isless(e1.reliability, e2.reliability)

function unwrap!(wrapped_image::AbstractArray{T, N};
                 range::Number=2*convert(T, pi),
                 wrap_around::NTuple{N, Bool}=tuple(fill(false, N)...),
                 seed::Int=-1) where {T, N}

    if seed != -1
        srand(seed)
    end

    range_T = convert(T, range)

    pixel_image = init_pixels(wrapped_image)
    calculate_reliability(pixel_image, wrap_around, range_T)
    edges = Edge{N}[]
    num_edges = _predict_num_edges(size(wrapped_image), wrap_around)
    sizehint!(edges, num_edges)
    for idx_dim=1:N
        populate_edges!(edges, pixel_image, idx_dim, wrap_around[idx_dim], range_T)
    end

    sort!(edges, alg=MergeSort)
    gather_pixels!(pixel_image, edges)
    unwrap_image!(wrapped_image, pixel_image, range_T)

    return wrapped_image
end

function _predict_num_edges(size_img, wrap_around)
    num_edges = 0
    for (size_dim, wrap_dim) in zip(size_img, wrap_around)
        num_edges += prod(size_img) * (size_dim-1) ÷ size_dim + wrap_dim * prod(size_img) ÷ size_dim
    end
    return num_edges
end

# function to broadcast
function init_pixels(wrapped_image::AbstractArray{T, N}) where {T, N}
    pixel_image = similar(wrapped_image, Pixel{T})
    @Threads.threads for i in eachindex(wrapped_image)
        @inbounds pixel_image[i] = Pixel(wrapped_image[i])
    end
    return pixel_image
end

function gather_pixels!(pixel_image, edges)
    for edge in edges
        @inbounds p1 = pixel_image[edge.pixel_1]
        @inbounds p2 = pixel_image[edge.pixel_2]
        merge_groups!(edge, p1, p2)
    end
end

function unwrap_image!(image, pixel_image, range)
    @Threads.threads for i in eachindex(image)
        @inbounds image[i] = range * pixel_image[i].periods + pixel_image[i].val
    end
end

function wrap_val(val, range)
    wrapped_val  = val
    wrapped_val += ifelse(val >  range/2, -range, zero(val))
    wrapped_val += ifelse(val < -range/2,  range, zero(val))
    return wrapped_val
end

function find_period(val_left, val_right, range)
    difference = val_left - val_right
    period  = 0
    period += ifelse(difference >  range/2, -1, 0)
    period += ifelse(difference < -range/2,  1, 0)
    return period
end

function merge_groups!(edge, pixel_1, pixel_2)
    if is_differentgroup(pixel_1, pixel_2)
        # pixel 2 is alone in group
        if is_pixelalone(pixel_2)
            merge_pixels!(pixel_1, pixel_2, -edge.periods)
        elseif is_pixelalone(pixel_1)
            merge_pixels!(pixel_2, pixel_1, edge.periods)
        else
            if is_bigger(pixel_1, pixel_2)
                merge_into_group!(pixel_1, pixel_2, -edge.periods)
            else
                merge_into_group!(pixel_2, pixel_1, edge.periods)
            end
        end
    end
end

@inline function is_differentgroup(p1::Pixel, p2::Pixel)
    return p1.head !== p2.head
end
@inline function is_pixelalone(pixel::Pixel)
    return pixel.head === pixel.last
end
@inline function is_bigger(p1::Pixel, p2::Pixel)
    return length(p1) ≥ length(p2)
end

function merge_pixels!(pixel_base::Pixel, pixel_target::Pixel, periods)
    pixel_base.head.groupsize += pixel_target.head.groupsize
    pixel_base.head.last.next = pixel_target.head
    pixel_base.head.last = pixel_target.head.last
    pixel_target.head = pixel_base.head
    pixel_target.periods = pixel_base.periods + periods
end

function merge_into_group!(pixel_base::Pixel, pixel_target::Pixel, periods)
    add_periods = pixel_base.periods + periods - pixel_target.periods
    pixel = pixel_target.head
    while pixel ≠ nothing
        # merge all pixels in pixel_target's group to pixel_base's group
        if pixel !== pixel_target
            pixel.periods += add_periods
            pixel.head = pixel_base.head
        end
        pixel = pixel.next
    end
    # assign pixel_target to pixel_base's group last
    merge_pixels!(pixel_base, pixel_target, periods)
end

function populate_edges!(edges, pixel_image::Array{T, N}, dim, connected, range) where {T, N}
    size_img       = collect(size(pixel_image))
    size_img[dim] -= 1
    idx_step       = fill(0, N)
    idx_step[dim] += 1
    idx_step_cart  = CartesianIndex{N}(idx_step...)
    idx_size       = CartesianIndex{N}(size_img...)
    for i in CartesianRange(idx_size)
        push!(edges, Edge{N}(pixel_image, i, i+idx_step_cart, range))
    end
    if connected
        idx_step        = fill(0, N)
        idx_step[dim]   = -size_img[dim]
        idx_step_cart   = CartesianIndex{N}(idx_step...)
        edge_begin      = fill(1, N)
        edge_begin[dim] = size(pixel_image)[dim]
        edge_begin_cart = CartesianIndex{N}(edge_begin...)
        for i in CartesianRange(edge_begin_cart, CartesianIndex(size(pixel_image)))
            push!(edges, Edge{N}(pixel_image, i, i+idx_step_cart, range))
        end
    end
end

function calculate_reliability(pixel_image::AbstractArray{T, N}, wrap_around, range) where {T, N}
    # get the shifted pixel indices in CartesinanIndex form
    # This gets all the nearest neighbors (CartesionIndex{N}() = one(CartesianIndex{N}))
    pixel_shifts = collect(CartesianRange(-CartesianIndex{N}(),
                                           CartesianIndex{N}()))
    size_img = size(pixel_image)
    # inner loop
    for i in CartesianRange((CartesianIndex{N}()+1), (CartesianIndex{N}(size_img)-1))
        @inbounds pixel_image[i].reliability = calculate_pixel_reliability(pixel_image, i, pixel_shifts, range)
    end

    for (idx_dim, connected) in enumerate(wrap_around)
        if connected
            # first border
            pixel_shifts_border = copy(pixel_shifts)
            for (idx_ps, ps) in enumerate(pixel_shifts_border)
                # if the pixel shift goes out of bounds, we make the shift wrap
                if ps[idx_dim] == 1
                    new_ps          = fill(0, N)
                    new_ps[idx_dim] = -size_img[idx_dim]+1
                    pixel_shifts_border[idx_ps] = CartesianIndex{N}(new_ps...)
                end
            end
            border_begin          = fill(2, N)
            border_begin[idx_dim] = size_img[idx_dim]
            border_begin          = CartesianIndex{N}(border_begin...)
            border_end            = collect(size_img)-1
            border_end[idx_dim]   = size_img[idx_dim]
            border_end = CartesianIndex{N}(border_end...)
            for i in CartesianRange(border_begin, border_end)
                @inbounds pixel_image[i].reliability = calculate_pixel_reliability(pixel_image, i, pixel_shifts_border, range)
            end
            # second border
            pixel_shifts_border = copy!(pixel_shifts_border, pixel_shifts)
            for (idx_ps, ps) in enumerate(pixel_shifts_border)
                # if the pixel shift goes out of bounds, we make the shift wrap, this time to the other side
                if ps[idx_dim] == -1
                    new_ps          = fill(0, N)
                    new_ps[idx_dim] = size_img[idx_dim]-1
                    pixel_shifts_border[idx_ps] = CartesianIndex{N}(new_ps...)
                end
            end
            border_begin          = fill(2, N)
            border_begin[idx_dim] = 1
            border_begin          = CartesianIndex{N}(border_begin...)
            border_end            = collect(size_img)-1
            border_end[idx_dim]   = 1
            border_end = CartesianIndex{N}(border_end...)
            for i in CartesianRange(border_begin, border_end)
                @inbounds pixel_image[i].reliability = calculate_pixel_reliability(pixel_image, i, pixel_shifts_border, range)
            end
        end
    end
end

function calculate_pixel_reliability(pixel_image::AbstractArray{T, N}, pixel_index, pixel_shifts, range) where {T, N}
    sum_val = zero(fieldtype(T, :val))
    for pixel_shift in pixel_shifts
        @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shift].val - pixel_image[pixel_index].val), range)^2
    end
    return sum_val
end

# specialized pixel reliability calculations for different N
function calculate_pixel_reliability(pixel_image::AbstractArray{T, 2}, pixel_index, pixel_shifts, range) where T
    @inbounds D1 = wrap_val(pixel_image[pixel_index+pixel_shifts[2]].val - pixel_image[pixel_index].val, range)
    @inbounds D2 = wrap_val(pixel_image[pixel_index+pixel_shifts[4]].val - pixel_image[pixel_index].val, range)
    @inbounds H  = wrap_val(pixel_image[pixel_index+pixel_shifts[6]].val - pixel_image[pixel_index].val, range)
    @inbounds V  = wrap_val(pixel_image[pixel_index+pixel_shifts[8]].val - pixel_image[pixel_index].val, range)
    return H*H + V*V + D1*D1 + D2*D2
end

function calculate_pixel_reliability(pixel_image::AbstractArray{T, 3}, pixel_index, pixel_shifts, range) where T
    sum_val = zero(fieldtype(T, :val))
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[1]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[2]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[3]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[4]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[5]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[6]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[7]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[8]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[9]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[10]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[11]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[12]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[13]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[15]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[17]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[18]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[19]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[20]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[21]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[22]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[23]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[24]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[25]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[26]].val - pixel_image[pixel_index].val, range))^2
    @inbounds sum_val += (wrap_val(pixel_image[pixel_index+pixel_shifts[27]].val - pixel_image[pixel_index].val, range))^2
    return sum_val
end

A = rand(10)
unwrap(A; range=10)
end
