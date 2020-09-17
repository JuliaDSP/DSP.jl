module Findpeaks

export findpeaks, PeakInfo

struct PeakInfo{T}
    height :: T
    left_diff :: T
    right_diff :: T
    prominence :: T
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
`findpeaks(y::Array{T},
x::Array{S}=collect(1:length(y))
;min_height::T=minimum(y), min_prom::T=minimum(y),
min_dist::S=0, threshold::T=0 ) where {T<:Real,S}`\n
Returns indices of local maxima (sorted from highest peaks to lowest)
in 1D array of real numbers. Similar to MATLAB's findpeaks().\n
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
"""
function findpeaks(
                   y :: AbstractVector{T},
                   x :: AbstractVector{S} = collect(1:length(y))
                   ;
                   kwargs...
                  ) where {T, S}

    if isempty(y)
        return (empty(x), empty(x, PeakInfo))
    end

    if length(x) != length(y)
        lx, ly = length(x), length(y)
        throw(ArgumentError("`x` and `y` need to have the same length: x($lx) y($ly)"))
    end

    min_height = get(kwargs, :min_height, minimum(y))
    min_prom   = get(kwargs, :min_prom,   zero(y[1]))
    min_dist   = get(kwargs, :min_dist,   zero(x[1]))
    threshold  = get(kwargs, :threshold,  zero(y[1]))

    peaks = 1:length(y) |> collect

    heights = y[peaks]
    inds2keep = heights .> min_height
    @apply_indices!(inds2keep, peaks, heights)

    inds2keep, ldiffs, rdiffs = in_threshold(peaks, y, threshold)
    @apply_indices!(inds2keep, peaks, ldiffs, rdiffs, heights)

    inds2keep, proms = with_prominence(peaks, y, min_prom)
    @apply_indices!(inds2keep, peaks, ldiffs, rdiffs, heights, proms)

    inds2keep = with_distance(peaks, x, y, min_dist)
    @apply_indices!(inds2keep, peaks, ldiffs, rdiffs, proms, heights)
    
    peak_info = [PeakInfo(heights[i], ldiffs[i], rdiffs[i], proms[i]) for i in 1:length(inds2keep)]
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

end # module
