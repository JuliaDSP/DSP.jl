module Findpeaks

export findpeaks, PeakInfo

struct PeakInfo{T}
    height :: T
    left_diff :: T
    right_diff :: T
    prominence :: T
    plateau_range :: UnitRange{Int}
end

"""
Overwrites arrays specified in `args` by `args[inds]`.
"""
macro apply_indices!(inds, args...)
    Expr(
         :block,
         [ :($(esc(args[i])) = $(esc(args[i]))[$(esc(inds))]) for i in 1:length(args) ]...)
end

"""
Macro specific to findpeaks function."
"""
macro on_empty_return(peaks, x)
    quote 
        if isempty($(esc(peaks)))
            return ((), ())
        end
    end
end

"""
`findpeaks(y::Array{T},
x::Array{S}=collect(1:length(y))
;min_height::T=minimum(y), min_prom::T=minimum(y),
min_dist::S=0, threshold::T=0 ) where {T<:Real,S}`\n

Finds peaks in 1D array. Similar to MATLAB's findpeaks().\n
Returns:\n
* tuple of peak positions found
* a tuple of `PeakInfo` structures with additional information about peaks.

*Arguments*:\n
`y` -- data\n
*Optional*:\n
`x` -- x-data\n
*Keyword*:\n
`min_height` -- minimal peak height\n
`min_prom` -- minimal peak prominence\n
`min_dist` -- minimal peak distance (keeping highest peaks)\n
`threshold` -- minimal difference (absolute value) between
 peak and neighboring points\n
`min_plateau_points` -- minimal number of points a plateau needs to have to be considered a peak
`max_plateau_points` -- maximal number of points a plateau can have to be considered a peak

The order of filtering is:
**height -> threshold -> prominence -> distance -> plateaus**

*Warning* `NaN` values are currently not supported and any `NaN`s should be filtered out by the user before applying this function.
"""
function findpeaks(
                   y :: AbstractVector{T},
                   x :: AbstractVector{S} = collect(1:length(y))
                   ;
                   kwargs...
                  ) where {T, S}

    @on_empty_return(y, x)

    if length(x) != length(y)
        lx, ly = length(x), length(y)
        throw(ArgumentError("`x` and `y` need to have the same length: x($lx) y($ly)"))
    end

    
    new_peaks, new_peak_info
    min_height = get(kwargs, :min_height, minimum(y))
    min_prom   = get(kwargs, :min_prom,   zero(y[1]))
    min_dist   = get(kwargs, :min_dist,   zero(x[1]))
    threshold  = get(kwargs, :threshold,  zero(y[1]))
    min_plateau_points  = get(kwargs, :min_plateau_points, zero(length(y)))
    max_plateau_points  = get(kwargs, :max_plateau_points, typemax(Int))

    peaks = 1:length(y) |> collect

    heights = y[peaks]
    inds2keep = heights .> min_height
    @apply_indices!(inds2keep, peaks, heights)
    @on_empty_return(peaks, x)

    inds2keep, ldiffs, rdiffs = in_threshold(peaks, y, threshold)
    @apply_indices!(inds2keep, peaks, ldiffs, rdiffs, heights)
    @on_empty_return(peaks, x)

    inds2keep, proms = with_prominence(peaks, y, min_prom)
    @apply_indices!(inds2keep, peaks, ldiffs, rdiffs, heights, proms)
    @on_empty_return(peaks, x)

    inds2keep = with_distance(peaks, x, y, min_dist)
    @apply_indices!(inds2keep, peaks, ldiffs, rdiffs, proms, heights)
    @on_empty_return(peaks, x)

    peak_info = [PeakInfo{T}(heights[i], ldiffs[i], rdiffs[i], proms[i], 1:1) for i in 1:length(inds2keep)]
    
    peaks, peak_info = group_plateaus(peaks, peak_info, min_plateau_points, max_plateau_points)
    (peaks = peaks, peak_info = peak_info)
end

function in_threshold(
                      peaks :: AbstractVector{Int},
                      y :: AbstractVector{T},
                      threshold :: T,
                     ) where {T <: Real}

    peak_len = length(peaks)

    inds2keep = 1:peak_len |> collect
    rdiffs = zeros(T, peak_len)
    ldiffs = zeros(T, peak_len)

    dy = diff(y)

    n_peaks2keep = 0
    for (j, peak_i) in enumerate(peaks)
        # resolve edges of the data
        if peak_i == 1
            # only check right difference
            rd = -dy[peak_i]
            if rd >= threshold
                n_peaks2keep += 1
                inds2keep[n_peaks2keep] = j
                ldiffs[j] = zero(threshold)
                rdiffs[j] = rd
            end
        elseif peak_i == length(y)
            # only check left difference
            ld = dy[peak_i - 1]
            if ld >= threshold
                n_peaks2keep += 1
                inds2keep[n_peaks2keep] = j
                ldiffs[j] = ld
                rdiffs[j] = zero(threshold)
            end
        else
            # check both sides
            ld = dy[peak_i - 1]
            rd = -dy[peak_i]
            if rd >= threshold && ld >= threshold
                n_peaks2keep += 1
                inds2keep[n_peaks2keep] = j
                ldiffs[j] = ld
                rdiffs[j] = rd
            end
        end
    end

    (inds2keep[1:n_peaks2keep], ldiffs, rdiffs)
end

function with_prominence(
                         peaks :: AbstractVector{Int},
                         y :: AbstractVector{T},
                         min_prom::T,
                        ) where {T <: Real}

    #minimal prominence refinement
    proms = prominence(y, peaks)
    proms .> min_prom, proms
end


function prominence(y :: AbstractVector{T},
                    peaks :: AbstractVector{Int}
                   ) where {T <: Real}
    yP = y[peaks]
    proms = zero(yP)

    for (i, p) in enumerate(peaks)
        lP, rP = 1, length(y)
        for j = (i-1):-1:1
            if yP[j] > yP[i]
                lP = peaks[j]
                break
            end
        end
        ml = minimum(y[lP:p])
        for j = (i+1):length(yP)
            if yP[j] > yP[i]
                rP = peaks[j]
                break
            end
        end
        mr = minimum(y[p:rP])
        ref = max(mr,ml)
        proms[i] = yP[i] - ref
    end

    proms
end

"""
Select only peaks that are further apart than `min_dist`
"""
function with_distance(
                       peaks :: AbstractVector{Int},
                       x :: AbstractVector{S},
                       y :: AbstractVector{T},
                       min_dist::S,
                      ) where {T <: Real, S}

    peaks2del = zeros(Bool, length(peaks))
    inds = sortperm(y[peaks], rev=true)
    sorted_peaks = copy(peaks)
    permute!(sorted_peaks, inds)
    for i = 1:length(sorted_peaks)
        for j = 1:(i-1)
            if abs(x[sorted_peaks[i]] - x[sorted_peaks[j]]) <= min_dist
                if !peaks2del[j]
                    peaks2del[i] = true
                end
            end
        end
    end

    inds[.!peaks2del]
end

function find_plateaus(peaks :: AbstractVector{Int})
    n_ranges = 0
    ranges = [1:1 for i in 1:length(peaks)]
    is_plateau = diff(peaks) .== 1
    constructed_range = -1:-1
    for (i, t) in enumerate(is_plateau)
        if t
            if constructed_range == -1:-1
                # start plateau range
                constructed_range = i:i
            end
        else
            if constructed_range != -1:-1
                # finish plateau range
                constructed_range = constructed_range.start : i
                n_ranges += 1
                ranges[n_ranges] = constructed_range

                constructed_range = -1:-1
            else
                # add single peak
                n_ranges += 1
                ranges[n_ranges] = i:i
            end
        end
    end
    # finish the last plateau
    if constructed_range != -1:-1
        constructed_range = constructed_range.start : length(peaks)
        n_ranges += 1
        ranges[n_ranges] = constructed_range
    else
        n_ranges += 1
        ranges[n_ranges] = length(peaks) : length(peaks)
    end
    ranges[1:n_ranges]
end


function merge_l_r_diff(l_edge :: PeakInfo, r_edge :: PeakInfo)
    (l_edge.left_diff, r_edge.right_diff)
end

function global_plateau_range(peaks :: AbstractVector{Int},
                              peak_info :: Vector{PeakInfo{T}},
                              r :: UnitRange{Int}) where T <: Real
    l_edge = peaks[r.start]
    r_edge = peaks[r.stop]
    plateau_range = l_edge : r_edge

    l_edge_info = peak_info[r.start]
    r_edge_info = peak_info[r.stop]
    left_diff, right_diff = merge_l_r_diff(l_edge_info, r_edge_info)
    merged_peak_info = PeakInfo(
                                l_edge_info.height,
                                left_diff,
                                right_diff,
                                l_edge_info.prominence,
                                plateau_range,
                               )

    (plateau_range, merged_peak_info)
end

function center_peak(plateau :: UnitRange{Int})
    center_peak = floor(Int, (plateau.start + plateau.stop)/2)
end

function group_plateaus(peaks :: AbstractVector{Int},
                        peak_info :: Vector{PeakInfo{T}},
                        min_plateau_points :: Int,
                        max_plateau_points :: Int) where T <: Real

    inds = sortperm(peaks)
    @apply_indices!(inds, peaks, peak_info)

    plateaus = find_plateaus(peaks)
    ranges = filter(x -> min_plateau_points <= length(x) <= max_plateau_points, plateaus)

    l = length(ranges)
    new_peaks = Vector{Int}(undef, l)
    new_peak_info = Vector{PeakInfo}(undef, l)
    for (i, r) in enumerate(ranges)
        plateau_range, merged_peak_info = global_plateau_range(peaks, peak_info, r)
        new_peaks[i] = center_peak(plateau_range)
        new_peak_info[i] = merged_peak_info
    end
    
    new_peaks, new_peak_info
end

end # module
