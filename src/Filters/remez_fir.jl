# linear phase FIR filter design that optimizes maximum error
# in the frequency domain

#=============================================

Translated from C code in scipy into Julia.
Tom Krauss, 2018.

SCIPY license:
Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2017 SciPy Developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of Enthought nor the names of the SciPy Developers
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

--- end of scipy license

The remez and associated functions were extracted from sigtoolsmodule.c,
with the following header comment:

SIGTOOLS module by Travis Oliphant

Copyright 2005 Travis Oliphant
Permission to use, copy, modify, and distribute this software without fee
is granted under the SciPy License.

C CODE BANNER

/********************************************************
 *
 *  Code taken from remez.c by Erik Kvaleberg which was
 *    converted from an original FORTRAN by:
 *
 * AUTHORS: JAMES H. MCCLELLAN
 *
 *         DEPARTMENT OF ELECTRICAL ENGINEERING AND COMPUTER SCIENCE
 *         MASSACHUSETTS INSTITUTE OF TECHNOLOGY
 *         CAMBRIDGE, MASS. 02139
 *
 *         THOMAS W. PARKS
 *         DEPARTMENT OF ELECTRICAL ENGINEERING
 *         RICE UNIVERSITY
 *         HOUSTON, TEXAS 77001
 *
 *         LAWRENCE R. RABINER
 *         BELL LABORATORIES
 *         MURRAY HILL, NEW JERSEY 07974
 *
 *
 *  Adaptation to C by
 *      egil kvaleberg
 *      husebybakken 14a
 *      0379 oslo, norway
 *  Email:
 *      egil@kvaleberg.no
 *  Web:
 *      http://www.kvaleberg.com/
 *
 *********************************************************/

=============================================#

# RemezFilterType:
#    Type I and II symmetric linear phase: neg==0   (filter_type==bandpass)
#    Type III and IV negative symmetric linear phase: neg==1   (filter_type==hilbert or differentiator)
@enum RemezFilterType filter_type_bandpass filter_type_differentiator filter_type_hilbert


# The strided product order limits roundoff in high-order interpolation.
function _remez_barycentric_weight(k, count, stride, nodes)
    product = 1.0
    node = nodes[k]
    @inbounds for offset in 1:stride, j in offset:stride:count
        j == k && continue
        product *= 2.0 * (node - nodes[j])
    end
    return inv(product)
end

_remez_callable(value::Real) = Returns(value)
_remez_callable(value) = value

function _remez_band_edges(band, Hz, lower, upper)
    low = Float64(clamp(band.first[1] / Hz, lower, upper))
    high = Float64(clamp(band.first[2] / Hz, lower, upper))
    return low, high
end

function _remez_grid(numtaps, band_defs, Hz, density, neg)
    odd = isodd(numtaps)
    ncosines = numtaps ÷ 2 + (odd && !neg)
    step = 0.5 / (density * ncosines)
    lower = neg ? step : 0.0
    upper = neg == odd ? 0.5 - step : 0.5

    # Count first, then fill directly; no normalized band or frequency arrays.
    ngrid = 0
    for band in band_defs
        low, high = _remez_band_edges(band, Hz, lower, upper)
        ngrid += max(length(low:step:high), 1)
    end
    grid = Vector{Float64}(undef, ngrid)
    desired = similar(grid)
    weight = similar(grid)

    index = 1
    for band in band_defs
        low, high = _remez_band_edges(band, Hz, lower, upper)
        response, importance = band.second isa Tuple{Any,Any} ?
                               band.second : (band.second, 1.0)
        index = _remez_grid_band!(grid, desired, weight, index, low:step:high,
                                  high, neg, odd, Hz,
                                  _remez_callable(response), _remez_callable(importance))
    end
    @assert index == ngrid + 1
    return grid, desired, weight
end

# Specialize on each band's callables, including heterogeneous band definitions.
function _remez_grid_band!(grid, desired, weight, index, frequencies, high, neg, odd, Hz,
                           response::D, importance::W) where {D,W}
    count = max(length(frequencies), 1)
    for i in 1:count
        # The final point is exactly the upper edge, including zero-width bands.
        frequency = i == count ? high : frequencies[i]
        change = neg ? sinpi(odd ? 2frequency : frequency) :
                       (odd ? 1.0 : cospi(frequency))
        grid[index] = cospi(2frequency)
        desired[index] = response(frequency * Hz) / change
        weight[index] = importance(frequency * Hz) * change
        index += 1
    end
    return index
end

# Barycentric evaluation away from the interpolation nodes. Callers use the
# stored values at nodes, so the hot loop needs no equality check or allocation.
@inline function _remez_eval(frequency, nodes, values, weights)
    denominator = 0.0
    numerator = 0.0
    # Independent divisions can run in parallel. Keep both accumulations in
    # their original order: reassociating these sums changes Remez convergence.
    j = 1
    @inbounds while j + 3 <= length(weights)
        c1 = weights[j] / (frequency - nodes[j])
        c2 = weights[j+1] / (frequency - nodes[j+1])
        c3 = weights[j+2] / (frequency - nodes[j+2])
        c4 = weights[j+3] / (frequency - nodes[j+3])
        denominator += c1
        numerator = fma(c1, values[j], numerator)
        denominator += c2
        numerator = fma(c2, values[j+1], numerator)
        denominator += c3
        numerator = fma(c3, values[j+2], numerator)
        denominator += c4
        numerator = fma(c4, values[j+3], numerator)
        j += 4
    end
    @inbounds for j in j:length(weights)
        coefficient = weights[j] / (frequency - nodes[j])
        denominator += coefficient
        numerator = muladd(coefficient, values[j], numerator)
    end
    return numerator / denominator
end

# Fit the alternating weighted error on the current extremal set.
function _remez_interpolate!(nodes, values, weights, extrema, grid, desired, importance)
    count = length(values)
    stride = (count - 2) ÷ 15 + 1
    for j in 1:count
        nodes[j] = grid[extrema[j]]
    end
    for j in 1:count
        weights[j] = _remez_barycentric_weight(j, count, stride, nodes)
    end

    numerator = denominator = 0.0
    sign = 1
    for j in 1:count
        index = extrema[j]
        numerator = muladd(weights[j], desired[index], numerator)
        denominator = muladd(sign, weights[j] / importance[index], denominator)
        sign = -sign
    end
    deviation = numerator / denominator
    first_sign = deviation > 0.0 ? -1 : 1
    ripple = -first_sign * deviation
    sign = first_sign
    for j in 1:count
        index = extrema[j]
        values[j] = desired[index] + sign * ripple / importance[index]
        sign = -sign
    end
    return first_sign, ripple
end

"""
    remez(numtaps::Integer, band_defs;
          Hz::Real=1.0,
          neg::Bool=false,
          maxiter::Integer=25,
          grid_density::Integer=16)

Calculate the minimax optimal filter using the Remez exchange algorithm [^McClellan1973a] [^McClellan1973b].

This is the simplified API that accepts just 2 required arguments (numtaps, band_defs).
For a scipy compatible version see the 3 arguments version (numtaps, bands, desired).

Calculate the filter-coefficients for the finite impulse response
(FIR) filter whose transfer function minimizes the maximum error
between the desired gain and the realized gain in the specified
frequency bands using the Remez exchange algorithm.

# Arguments
- `numtaps::Integer`: The desired number of taps in the filter.
    The number of taps is the number of terms in the filter, or the filter
    order plus one.
- `band_defs`: A sequence of band definitions.
    This sequence defines the bands. Each entry is a pair. The pair's
    first item is a tuple of band edges (low, high). The pair's second item
    defines the desired response and weight in that band. The weight is optional
    and defaults to 1.0. Both the desired response and weight may be either scalars
    or functions. If a function, the function should accept a real frequency and
    return the real desired response or real weight. Examples:
    + LPF with unity weights. `[(0, 0.475) => 1, (0.5, 1.0) => 0]`
    + LPF with weight of 2 in the stop band. `[(0, 0.475) => (1, 1), (0.5, 1.0) => (0, 2)]`
    + BPF with unity weights. `[(0, 0.375) => 0, (0.4, 0.5) => 1, (0.525, 1.0) => 0]`
    + Hilbert transformer. `[(0.1, 0.95) => 1]; neg=true`
    + Differentiator. `[(0.01, 0.99) => (f -> f/2, f -> 1/f)]; neg=true`
- `Hz::Real`: The sampling frequency in Hz. Default is 1.
- `neg::Bool`: Whether the filter has negative symmetry or not. Default is false.
    If false, the filter is even-symmetric. If true, the filter is odd-symmetric.
    neg=true means that h[n]=-h[end+1-n]; neg=false means that h[n]=h[end+1-n].
- `maxiter::Integer`: (optional)
    Maximum number of iterations of the algorithm. Default is 25.
- `grid_density:Integer`: (optional)
    Grid density. The dense grid used in `remez` is of size
    `(numtaps + 1) * grid_density`. Default is 16.

# Returns
- `h::Array{Float64,1}`: A rank-1 array containing the coefficients of the optimal
    (in a minimax sense) filter.

[^McClellan1973a]: J. H. McClellan and T. W. Parks,
    A unified approach to the design of optimum FIR linear phase digital filters,
    IEEE Trans. Circuit Theory, vol. CT-20, pp. 697-701, 1973.

[^McClellan1973b]: J. H. McClellan, T. W. Parks and L. R. Rabiner,
    A Computer Program for Designing Optimum FIR Linear Phase Digital Filters,
    IEEE Trans. Audio Electroacoust., vol. AU-21, pp. 506-525, 1973.

# Examples
Construct a length 35 filter with a passband at 0.15-0.4 Hz
(desired response of 1), and stop bands at 0-0.1 Hz and 0.45-0.5 Hz
(desired response of 0). Note: the behavior in the frequency ranges between
those bands - the transition bands - is unspecified.

```jldoctest
julia> bpass = remez(35, [(0, 0.1)=>0, (0.15, 0.4)=>1, (0.45, 0.5)=>0]);
```

You can trade-off maximum error achieved for transition bandwidth.
The wider the transition bands, the lower the maximum error in the
bands specified. Here is a bandpass filter with the same passband, but
wider transition bands.

```jldoctest
julia> bpass2 = remez(35, [(0, 0.08)=>0, (0.15, 0.4)=>1, (0.47, 0.5)=>0]);
```

Here we compute the frequency responses and plot them in dB.

```julia-repl
julia> using PyPlot
julia> b = PolynomialRatio(bpass, [1.0])
julia> b2 = PolynomialRatio(bpass2, [1.0])
julia> f = range(0, stop=0.5, length=1000)
julia> plot(f, 20*log10.(abs.(freqresp(b,f,1.0))))
julia> plot(f, 20*log10.(abs.(freqresp(b2,f,1.0))))
julia> grid()
```

# Examples from the unittests - standard (even) symmetry.

Length 151 LPF (Low Pass Filter).
```jldoctest
julia> h = remez(151, [(0, 0.475) => 1, (0.5, 1.0) => 0]; Hz=2.0);
```

Length 152 LPF. Non-default "weight" input.
```jldoctest
julia> h = remez(152, [(0, 0.475) => (1, 1), (0.5, 1.0) => (0, 2)]; Hz=2.0);
```

Length 51 HPF (High Pass Filter).
```jldoctest
julia> h = remez(51, [(0, 0.75) => 0, (0.8, 1.0) => 1]; Hz=2.0);
```

Length 180 BPF (Band Pass Filter).
```jldoctest
julia> h = remez(180, [(0, 0.375) => 0, (0.4, 0.5) => 1, (0.525, 1.0) => 0]; Hz=2.0, maxiter=30);
```

# Examples from the unittests - Odd-symmetric filters - hilbert and differentiators type.
Even length - has a much better approximation since the response is not constrained to 0 at
the nyquist frequency.  Length 20 Hilbert transformer.
```jldoctest
julia> h = remez(20, [(0.1, 0.95) => 1]; neg=true, Hz=2.0);
```

Length 21 Hilbert transformer.
```jldoctest
julia> h = remez(21, [(0.1, 0.95) => 1]; neg=true, Hz=2.0);
```

Length 200 differentiator.
```jldoctest
julia> h = remez(200, [(0.01, 0.99) => (f -> f/2, f -> 1/f)]; neg=true, Hz=2.0);
```

Length 201 differentiator.
```jldoctest
julia> h = remez(201, [(0.05, 0.95) => (f -> f/2, f -> 1/f)]; neg=true, Hz=2.0);
```

Inverse sinc filter - custom response function
```julia-repl
julia> L = 64; Fs = 4800*L;
julia> passband_response_function = f -> (f==0) ? 1.0 : abs.((π*f/4800) ./ sin.(π*f/4800));
julia> h = remez(201, [(    0.0, 2880.0) => (passband_response_function, 1.0),
                (10000.0,   Fs/2) => (0.0, 100.0)]; Hz=Fs);
```
"""
function remez(numtaps::Integer, band_defs;
                Hz::Real=1.0,
                neg::Bool=false,
                maxiter::Integer=25,
                grid_density::Integer=16)
    all(b.first[1] <= b.first[2] for b in band_defs) ||
        throw(ArgumentError("lower band edge higher then upper band edge"))
    all(band_defs[i].first[2] < band_defs[i+1].first[1] for i in 1:length(band_defs)-1) ||
        throw(ArgumentError("band edges is not monotonically increasing"))
    (0 <= band_defs[1].first[1]) && (band_defs[end].first[2] <= 0.5*Hz) ||
        throw(ArgumentError("band edges must be between 0 and `Hz`/2"))

    grid, desired, weight = _remez_grid(numtaps, band_defs, Hz, grid_density, neg)
    full_grid = band_defs[1].first[1] == 0.0 && band_defs[end].first[2] == 0.5Hz
    return _remez_solve(numtaps, neg, maxiter, full_grid, grid, desired, weight)
end

function _remez_solve(numtaps, neg, maxiter, full_grid, grid, desired, weight)
    ncosines = numtaps ÷ 2 + (isodd(numtaps) && !neg)
    count = ncosines + 1
    extrema = Vector{Int}(undef, count + 1)
    nodes = zeros(count + 1) # Extra node is a sentinel during coefficient recovery.
    values = zeros(count)
    weights = zeros(count)
    for j in 1:count
        extrema[j] = (j - 1) * (length(grid) - 1) ÷ ncosines + 1
    end
    extrema[end] = length(grid) + 1

    previous_ripple = -1.0
    iteration = 0
    while true
        iteration += 1
        if iteration > maxiter
            @warn("remez() iteration count exceeds maxiter = $maxiter, filter is not converged; try increasing maxiter")
            break
        end

        first_sign, ripple = _remez_interpolate!(nodes, values, weights, extrema,
                                                 grid, desired, weight)
        if ripple <= previous_ripple
            throw(ErrorException("remez() - failure to converge at iteration $iteration, try reducing transition band width"))
        end
        previous_ripple = ripple
        _remez_exchange!(extrema, grid, desired, weight, nodes, values, weights,
                         first_sign, ripple) || break
    end

    return _remez_coefficients(numtaps, neg, full_grid, grid, nodes, values, weights)
end

# Search a monotone run of signed error, stopping before the exclusive bound.
# The direction is a value parameter so the two hot loops specialize separately.
@inline function _remez_peak(error, index, bound, amplitude, ::Val{step}) where {step}
    while index + step != bound
        candidate = error(index + step)
        candidate <= amplitude && break
        index += step
        amplitude = candidate
    end
    return index, amplitude
end

# Exchange extrema in place. Searching each interval in the same order preserves
# the original algorithm's tie breaking and its behavior at the iteration limit.
function _remez_exchange!(extrema, grid, desired, weight, nodes, values, baryweights,
                          first_sign, ripple)
    count = length(values)
    gridlength = length(grid)
    oldfirst, oldlast = extrema[1], extrema[count]
    lower = 0
    changed = false
    first_amplitude = last_amplitude = ripple
    sign = first_sign

    for j in 1:count
        current, upper = extrema[j], extrema[j + 1]
        error = let sign = sign
            i -> sign * ((_remez_eval(grid[i], nodes, values, baryweights) -
                         desired[i]) * weight[i])
        end
        amplitude = ripple
        candidate = current + 1

        value = candidate < upper ? error(candidate) : amplitude
        if value > amplitude
            candidate, amplitude = _remez_peak(error, candidate, upper,
                                               value, Val(1))
            lower = candidate
            changed = true
        else
            candidate = current - 1
            while candidate > lower
                value = error(candidate)
                if value > amplitude
                    candidate, amplitude = _remez_peak(error, candidate, lower,
                                                       value, Val(-1))
                    break
                end
                changed && break
                candidate -= 1
            end

            if candidate > lower && amplitude > ripple
                # Keep the old extremum as the lower bound when moving left.
                lower = current
                changed = true
            elseif candidate <= lower && changed
                candidate = current
                lower = current
            elseif !changed
                candidate = current + 2
                while candidate < upper
                    value = error(candidate)
                    value > amplitude && break
                    candidate += 1
                end
                if candidate < upper
                    candidate, amplitude = _remez_peak(error, candidate, upper,
                                                       value, Val(1))
                    changed = true
                else
                    candidate = current
                end
                lower = candidate
            else
                candidate = current
                lower = current
            end
        end

        extrema[j] = candidate
        j == 1 && (first_amplitude = amplitude)
        last_amplitude = amplitude
        sign = -sign
    end

    # An extra extremum at either end may displace the opposite endpoint.
    # The small margin avoids cycling on nearly equal endpoint errors.
    left_bound = min(oldfirst, extrema[1])
    right_bound = max(oldlast, extrema[count])
    left_error = i -> -first_sign * ((_remez_eval(grid[i], nodes, values, baryweights) -
                                     desired[i]) * weight[i])
    left = 1
    threshold = last_amplitude * 1.00001
    amplitude = threshold
    while left < left_bound
        amplitude = left_error(left)
        amplitude > threshold && break
        left += 1
    end
    have_left = left < left_bound
    if have_left
        left, amplitude = _remez_peak(left_error, left, left_bound,
                                      amplitude, Val(1))
        first_amplitude = max(first_amplitude, amplitude)
    end

    right_error = let sign = sign
        i -> sign * ((_remez_eval(grid[i], nodes, values, baryweights) -
                     desired[i]) * weight[i])
    end
    right = gridlength
    threshold = first_amplitude * 1.00001
    while right > right_bound
        amplitude = right_error(right)
        amplitude > threshold && break
        right -= 1
    end
    if right > right_bound
        right, _ = _remez_peak(right_error, right, right_bound,
                               amplitude, Val(-1))
        for j in 1:count-1
            extrema[j] = extrema[j + 1]
        end
        extrema[count] = right
        return true
    elseif have_left
        for j in count:-1:2
            extrema[j] = extrema[j - 1]
        end
        extrema[1] = left
        return true
    end
    return changed
end

function _remez_coefficients(numtaps, neg, full_grid, grid, nodes, values, weights)
    ncosines = length(values) - 1
    count, workspace_length = ncosines + 1, ncosines + 2
    gridlength = length(grid)
    # Sample the interpolant, then recover its cosine coefficients by an IDFT.

    samples = Vector{Float64}(undef, ncosines) # frequency response on evenly spaced grid
    coefficients = zeros(Float64, workspace_length)   # return vector

    node_tolerance = 1.0e-06
    nodes[workspace_length] = -2.0
    step = 1 / (2*ncosines - 1)
    node_index = 1

    # Partial grids use an affine change of the Chebyshev variable.
    full_grid |= ncosines <= 3
    if !full_grid
        scale    = 2.0/(grid[1]-grid[gridlength])
        shift    = -(grid[1]+grid[gridlength])/(grid[1]-grid[gridlength])
    end

    # Reuse values near interpolation nodes to avoid division by zero.
    for j = 1:ncosines
        frequency = cospi(2 * (j - 1) * step)
        if !full_grid
            frequency = (frequency-shift)/scale
        end
        if (node_index > 1)
            node_index = node_index-1
        end
        while nodes[node_index]-frequency >= node_tolerance
            node_index += 1
        end
        if frequency-nodes[node_index] < node_tolerance
            samples[j] = values[node_index]
        else
            samples[j] = _remez_eval(frequency, nodes, values, weights)
        end
    end

    degree = ncosines - 1

    for j in 1:ncosines
        total = 0.0
        for k in 1:degree
            total = muladd(samples[k+1], cospi(2 * (j-1) * step * k), total)
        end
        coefficients[j] = 2total + samples[1]
    end

    for j in 2:ncosines
        coefficients[j] *= 2.0 * step
    end
    coefficients[1] *= step

    if !full_grid
        # Interpolation is finished; reuse its workspace for basis conversion.
        p = weights
        q = values
        p[1] = muladd(2coefficients[ncosines], shift, coefficients[degree])
        p[2] = 2.0*scale*coefficients[ncosines]
        q[1] = coefficients[ncosines-2]-coefficients[ncosines]
        for j in 2:degree
            if j >= degree
                scale *= 0.5
                shift *= 0.5
            end
            p[j+1] = 0.0
            for k in 1:j
                samples[k] = p[k]
                p[k] = 2.0 * shift * samples[k]
            end
            p[2] = muladd(samples[1], 2scale, p[2])
            for k in 1:j-1
                p[k] += muladd(scale, samples[k+1], q[k])
            end
            for k in 3:j+1
                p[k] = muladd(scale, samples[k-1], p[k])
            end

            if j != degree
                for k in 1:j
                    q[k] = -samples[k]
                end
                q[1] += coefficients[ncosines - 1 - j]
            end
        end
        for j in 1:ncosines
            coefficients[j] = p[j]
        end
    end

    return _remez_impulse_response(numtaps, neg, coefficients, ncosines)
end

# Undo the linear-phase transformation for FIR types I, II, III and IV.
function _remez_impulse_response(numtaps, neg, coefficients, ncosines)
    odd = isodd(numtaps)
    count = ncosines + 1
    degree = ncosines - 1
    h = Vector{Float64}(undef, numtaps)
    if !neg
        if odd
            for j in 1:degree
                h[j] = 0.5 * coefficients[count-j]
            end
            h[ncosines] = coefficients[1]
        else
            h[1] = 0.25 * coefficients[ncosines]
            for j in 2:degree
                h[j] = 0.25 * (coefficients[count-j] + coefficients[ncosines+2-j])
            end
            h[ncosines] = muladd(0.5, coefficients[1], 0.25 * coefficients[2])
        end
    else
        if odd
            for j in 1:degree
                h[j] = 0.25 * (coefficients[count-j] - coefficients[ncosines+3-j])
            end
            h[ncosines] = muladd(0.5, coefficients[1], -0.25 * coefficients[3])
        else
            h[1] = 0.25 * coefficients[ncosines]
            for j in 2:degree
                h[j] = 0.25 * (coefficients[count-j] - coefficients[ncosines+2-j])
            end
            h[ncosines] = muladd(0.5, coefficients[1], -0.25 * coefficients[2])
        end
    end

    for j in 1:ncosines
        k = numtaps + 1 - j
        if !neg
           h[k] = h[j]
        else
           h[k] = -h[j]
        end
    end
    if neg && odd
        h[count] = 0.0
    end

    return h
end


"""
    remez(numtaps::Integer,
          bands::Vector,
          desired::Vector;
          weight::Vector=[],
          Hz::Real=1.0,
          filter_type::RemezFilterType=filter_type_bandpass,
          maxiter::Integer=25,
          grid_density::Integer=16)

This is the scipy compatible version that requires 3 arguments (numtaps, bands, desired).
For a simplified API, see the 2 argument version (numtaps, band_defs). The filters
designed are equivalent, the inputs are just specified in a different way.
Below the arguments and examples are described that differ from the simplified
API version.

# Arguments
- `bands::Vector`: A monotonic sequence containing the band edges in Hz.
    All elements must be non-negative and less than half the sampling
    frequency as given by `Hz`.
- `desired::Vector`:A sequence half the size of bands containing the desired
    gain in each of the specified bands.
- `weight::Vector`: (optional)
    A relative weighting to give to each band region. The length of
    `weight` has to be half the length of `bands`.
- `filter_type::RemezFilterType`: Default is `filter_type_bandpass`.
    The type of filter:
    +  `filter_type_bandpass` : flat response in bands. This is the default.
    +  `filter_type_differentiator` : frequency proportional response in bands.
        Odd symetric as in `filter_type_hilbert` case, but with a linear sloping
        desired response.
    +  `filter_type_hilbert` : filter with odd symmetry, that is, type III
                  (for even order) or type IV (for odd order)
                  linear phase filters.

# Examples
Compare the examples with the simplified API and the Scipy API.
Each of the following blocks first designs a filter using the
simplified (recommended) API, and then designs the same filter
using the Scipy-compatible API.

```jldoctest
julia> bpass = remez(35, [(0, 0.1)=>0, (0.15, 0.4)=>1, (0.45, 0.5)=>0]);

julia> bpass = remez(35, [0, 0.1, 0.15, 0.4, 0.45, 0.5], [0, 1, 0]);

```

```jldoctest
julia> bpass2 = remez(35, [(0, 0.08)=>0, (0.15, 0.4)=>1, (0.47, 0.5)=>0]);

julia> bpass2 = remez(35, [0, 0.08, 0.15, 0.4, 0.47, 0.5], [0, 1, 0]);

```

```jldoctest
julia> h = remez(20, [(0.1, 0.95) => 1]; neg=true, Hz=2.0);

julia> h = remez(20, [0.1, 0.95], [1]; filter_type=filter_type_hilbert, Hz=2.0);

```

```jldoctest
julia> h = remez(200, [(0.01, 0.99) => (f -> f/2, f -> 1/f)]; neg=true, Hz=2.0);

julia> h = remez(200, [0.01, 0.99], [1]; filter_type=filter_type_differentiator, Hz=2.0);

```
"""
function remez(numtaps::Integer, bands::Vector, desired::Vector;
               weight::Vector=fill(1.0, length(desired)),
               Hz::Real=1.0,
               filter_type::RemezFilterType=filter_type_bandpass,
               kwargs...)
    issorted(bands) || throw(ArgumentError("`bands` is not monotonically increasing"))
    length(bands) == 2length(desired) ||
        throw(ArgumentError("`desired` must be half the length of `bands`."))
    length(bands) == 2length(weight) ||
        throw(ArgumentError("`weight` must be half the length of `bands`."))
    band_ranges = [(bands[i], bands[i+1]) for i in 1:2:length(bands)]
    if filter_type == filter_type_differentiator
        eff = [f -> d*f/Hz for d in desired]
        wate = [d > 0.0001 ? f -> w/(f/Hz) : f -> w for (w, d) in zip(weight, desired)]
    else
        eff = desired
        wate = weight
    end
    band_defs = [r => (d, w) for (r, d, w) in zip(band_ranges, eff, wate)]
    neg = filter_type in (filter_type_hilbert, filter_type_differentiator)
    return remez(numtaps, band_defs; Hz, neg, kwargs...)
end
