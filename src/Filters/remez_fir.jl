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



"""
    lagrange_interp(k::Integer, n::Integer, m::Integer, x::AbstractVector)

CALCULATE THE LAGRANGE INTERPOLATION COEFFICIENTS
"""
function lagrange_interp(k::Integer, n::Integer, m::Integer, x::AbstractVector)
    retval = 1.0
    q = x[k]
    for l = 1 : m
        for j = l : m : n
            if j != k
                retval *= 2.0 * (q - x[j])
            end
        end
    end
    1.0 / retval
end


"""
    build_grid(numtaps, band_defs, Hz, grid_density, neg)

Returns `grid`, `des`, and `wt` arrays
"""
function build_grid(nfilt, band_defs, Hz, grid_density, neg)
    nodd = isodd(nfilt)
    nfcns = nfilt ÷ 2
    if nodd && !neg
        nfcns = nfcns + 1
    end

    #
    # SET UP THE DENSE GRID. THE NUMBER OF POINTS IN THE GRID
    # IS (FILTER LENGTH + 1)*GRID DENSITY/2
    #
    delf = 0.5 / (grid_density * nfcns)

    flimlow = neg ? delf : 0.0
    flimhigh = (neg == nodd) ? 0.5 - delf : 0.5
    function normalize_banddef_entry(b)
        # normalize and clamp band-edges
        fl = convert(Float64, clamp(b.first[1] / Hz, flimlow, flimhigh))
        fu = convert(Float64, clamp(b.first[2] / Hz, flimlow, flimhigh))
        # make sure desired and weight are functions
        if !isa(b.second, Tuple{Any, Any})
            desired = b.second
            weight = 1.0
        else
            desired = b.second[1]
            weight = b.second[2]
        end
        if isa(desired, Real)
            let d = desired
                desired = _ -> d
            end
        end
        if isa(weight, Real)
            let w = weight
                weight = _ -> w
            end
        end
        return Pair{Tuple{Float64,Float64},Tuple{Any,Any}}((fl, fu), (desired, weight))
    end
    normalized_band_defs = normalize_banddef_entry.(band_defs)

    local ngrid
    # work around JuliaLang/julia#15276
    let delf=delf
        ngrid = sum(max(length(band_def.first[1]:delf:band_def.first[2]), 1) for band_def in normalized_band_defs)#::Int
    end

    grid = zeros(Float64, ngrid)  # the array of frequencies, between 0 and 0.5
    des = zeros(Float64, ngrid)   # the desired function on the grid
    wt = zeros(Float64, ngrid)    # array of weights

    j = 1

    #
    # CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT
    # FUNCTION ON THE GRID
    #
    # and
    #
    # SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT
    # TO THE ORIGINAL PROBLEM
    #
    for band_def in normalized_band_defs
        flow = band_def.first[1]
        fup = band_def.first[2]
        # outline inner loop to have it type-stable (band_def.second is Tuple{Any,Any})
        j = _buildgrid!(grid, des, wt, j, [(flow:delf:fup)[1:end-1]; fup], neg, nodd, Hz,
                        band_def.second)
    end
    @assert ngrid == j - 1

    return grid, des, wt
end

function _buildgrid!(grid, des, wt, j, fs, neg, nodd, Hz, des_wt)
    for f in fs
        change = neg ? (nodd ? sinpi(2f) : sinpi(f)) : (nodd ? 1.0 : cospi(f))
        grid[j] = cospi(2f)
        des[j] = des_wt[1](f * Hz) / change
        wt[j] = des_wt[2](f * Hz) * change
        j += 1
    end
    return j
end

"""
    function freq_eval(xf, x::AbstractVector, y::AbstractVector, ad::AbstractVector)
-----------------------------------------------------------------------
 FUNCTION: freq_eval (gee)
  FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE
  LAGRANGE INTERPOLATION FORMULA IN THE BARYCENTRIC FORM

-----------------------------------------------------------------------
"""
function freq_eval(xf, x::AbstractVector, y::AbstractVector, ad::AbstractVector)
    d = 0.0
    p = 0.0

    for j in eachindex(ad)
        c = ad[j] / (xf - x[j])
        d += c
        p = muladd(c, y[j], p)
    end

    p/d
end

function initialize_y(dev::Float64, nz::Integer, iext::AbstractArray, des::AbstractArray, wt::AbstractArray, y::AbstractArray)
    nu = 1
    if dev > 0.0
        nu = -1
    end
    dev = -nu * dev
    k = nu
    for j = 1 : nz
        l = iext[j]
        y[j] = des[l] + k * dev / wt[l]
        k = -k
    end
    nu, dev
end

#=========
Banner from C code

-----------------------------------------------------------------------
 SUBROUTINE: remez
  THIS SUBROUTINE IMPLEMENTS THE REMEZ EXCHANGE ALGORITHM
  FOR THE WEIGHTED CHEBYSHEV APPROXIMATION OF A CONTINUOUS
  FUNCTION WITH A SUM OF COSINES.  INPUTS TO THE SUBROUTINE
  ARE A DENSE GRID WHICH REPLACES THE FREQUENCY AXIS, THE
  DESIRED FUNCTION ON THIS GRID, THE WEIGHT FUNCTION ON THE
  GRID, THE NUMBER OF COSINES, AND AN INITIAL GUESS OF THE
  EXTREMAL FREQUENCIES.  THE PROGRAM MINIMIZES THE CHEBYSHEV
  ERROR BY DETERMINING THE BSMINEST LOCATION OF THE EXTREMAL
  FREQUENCIES (POINTS OF MAXIMUM ERROR) AND THEN CALCULATES
  THE COEFFICIENTS OF THE BEST APPROXIMATION.
-----------------------------------------------------------------------
=========#
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
- `bands_defs`: A sequence of band definitions.
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
julia> b = DSP.Filters.PolynomialRatio(bpass, [1.0])
julia> b2 = DSP.Filters.PolynomialRatio(bpass2, [1.0])
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

    grid, des, wt = build_grid(numtaps, band_defs, Hz, grid_density, neg)

    nfilt = numtaps
    ngrid = length(grid)

    nodd = isodd(nfilt)   # boolean: "nodd" means filter length is odd
    nfcns = numtaps ÷ 2   # integer divide
    if nodd && !neg
        nfcns = nfcns + 1
    end

    nz  = nfcns+1
    nzz = nfcns+2
    iext = zeros(Int64, nzz)   # indices of extremals
    x = zeros(Float64, nzz)
    y = zeros(Float64, nz)

    for j = 1:nz
        iext[j] = (j-1)*(ngrid-1) ÷ nfcns + 1
    end

    dev = 0.0     # deviation from the desired function,
                  # that is, the amount of "ripple" on the extremal set
    devl = -1.0   # deviation on last iteration
    niter = 0
    ad = zeros(Float64, nz)

    jet = ((nfcns-1) ÷ 15) + 1

    while true

        #
        # Start next iteration
        #
    #   @label L100
        iext[nzz] = ngrid + 1
        niter += 1
        if niter > maxiter
            @warn("remez() iteration count exceeds maxiter = $maxiter, filter is not converged; try increasing maxiter")
            # the filter is returned in its current, unconverged state.
            break
        end

        for j = 1:nz
            x[j] = grid[iext[j]]
        end

        for j = 1 : nz
            ad[j] = lagrange_interp(j, nz, jet, x)
        end

        dnum = 0.0
        dden = 0.0
        k = 1
        for j = 1 : nz
            l = iext[j]
            dnum = muladd(ad[j], des[l], dnum)
            dden = muladd(k, ad[j] / wt[l], dden)
            k = -k
        end
        dev = dnum / dden

        fill!(y, 0.0)
        nu, dev = initialize_y(dev, nz, iext, des, wt, y)

        if dev <= devl
            # finished
            throw(ErrorException("remez() - failure to converge at iteration $niter, try reducing transition band width"))
        end
        devl = dev

        #
        # SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST APPROXIMATION
        #

        # Between here and L370, the extremal index set is updated in a loop
        # roughly over the index "j" - although the logic is complicated as
        # the extremal set may grow or shrink in an iteration.
        # j - the index of the current extremal being updated
        # nz - the number of cosines in the approximation (including the constant term).
        #      nz = nfcns + 1 where nfcns = nfilt / 2, and
        #      nfilt is the filter length or number of taps.
        #      For example, for a length 15 filter, nfcns = 7 and nz = 8.
        # jchgne - the number of extremal indices that changed this iteration
        jchnge = 0
        k1 = iext[1]
        knz = iext[nz]
        klow = 0
        nut = -nu
        j = 1

        local comp

      @label L200
        j == nzz && (ynz = comp)   # equivalent to "if (j == nzz) ynz = comp; end"
        j >= nzz && @goto L300
        kup = iext[j+1]
        l = iext[j]+1
        nut = -nut
        j == 2 && (y1 = comp)
        comp = dev
        l >= kup && @goto L220
        err = (freq_eval(grid[l], x, y, ad) - des[l]) * wt[l]
        nut*err <= comp && @goto L220
        comp = nut * err
      @label L210
        l += 1; l >= kup && @goto L215
        err = (freq_eval(grid[l], x, y, ad) - des[l]) * wt[l]
        nut*err <= comp && @goto L215
        comp = nut * err
        @goto L210

      @label L215
        iext[j] = l - 1; j += 1
        klow = l - 1
        jchnge += 1
        @goto L200

      @label L220
        l -= 1
      @label L225
        l -= 1; l <= klow && @goto L250
        err = (freq_eval(grid[l], x, y, ad) - des[l]) * wt[l]
        nut*err > comp && @goto L230
        jchnge <= 0 && @goto L225
        @goto L260

      @label L230
        comp = nut * err
      @label L235
        l -= 1; l <= klow && @goto L240
        err = (freq_eval(grid[l], x, y, ad) - des[l]) * wt[l]
        nut*err <= comp && @goto L240
        comp = nut * err
        @goto L235
      @label L240
        klow = iext[j]
        iext[j] = l+1
        j += 1
        jchnge += 1
        @goto L200

      @label L250
        l = iext[j]+1
        jchnge > 0 && @goto L215

      @label L255
        l += 1; l >= kup && @goto L260
        err = (freq_eval(grid[l], x, y, ad) - des[l]) * wt[l]
        nut*err <= comp && @goto L255
        comp = nut * err

        @goto L210
      @label L260
        klow = iext[j]; j += 1
        @goto L200

      @label L300
        j > nzz && @goto L320
        k1 > iext[1] && (k1 = iext[1])
        knz < iext[nz] && (knz = iext[nz])
        nut1 = nut
        nut = -nu
        l = 0
        kup = k1
        comp = ynz*(1.00001)
        luck = 1
      @label L310
        l += 1; l >= kup && @goto L315
        err = (freq_eval(grid[l], x, y, ad) - des[l]) * wt[l]
        nut*err <= comp && @goto L310
        comp =  nut * err
        j = nzz
        @goto L210

      @label L315
        luck = 6
        @goto L325

      @label L320
        luck > 9 && @goto L350
        comp > y1 && (y1 = comp)
        k1 = iext[nzz]
      @label L325
        l = ngrid+1
        klow = knz
        nut = -nut1
        comp = y1*(1.00001)
      @label L330
        l -= 1; l <= klow && @goto L340
        err = (freq_eval(grid[l], x, y, ad) - des[l]) * wt[l]
        nut*err <= comp && @goto L330
        j = nzz
        comp =  nut * err
        luck = luck + 10
        @goto L235
      @label L340
        luck == 6 && @goto L370
        for j = 1 : nfcns
            iext[nzz-j] = iext[nz-j]
        end
        iext[1] = k1
        continue    # @goto L100
      @label L350
        for j = 1:nz
            iext[j] = iext[j+1]
        end

        continue    # @goto L100
      @label L370


        if jchnge <= 0  # we are done if none of the extremal indices changed
            break
        end
    end  # while

    #
    #    CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION
    #    USING THE INVERSE DISCRETE FOURIER TRANSFORM
    #

    a = zeros(Float64, nfcns)   # frequency response on evenly spaced grid
    p = zeros(Float64, nfcns)
    q = zeros(Float64, nfcns-2)
    alpha = zeros(Float64, nzz)   # return vector

    fsh = 1.0e-06
    x[nzz] = -2.0
    delf = 1 / (2*nfcns - 1)
    l = 1

    # Boolean for "kkk" in C code.
    full_grid = (band_defs[1].first[1] == 0.0 && band_defs[end].first[2] == 0.5*Hz) || (nfcns <= 3)
    if !full_grid
        aa    = 2.0/(grid[1]-grid[ngrid])
        bb    = -(grid[1]+grid[ngrid])/(grid[1]-grid[ngrid])
    end

    # Fill in "a" array with the frequency response on an evenly
    # spaced set of frequencies. Care is taken to use "y[l]" -
    # an already computed response on one of the extremals "l" -
    # if the extremal is equal to the frequency ft. If no y[l]
    # matches, a[j] is computed using freq_eval.
    for j = 1 : nfcns
        xt = cospi(2 * (j - 1) * delf)
        if !full_grid
            xt = (xt-bb)/aa
        end
        if (l > 1)
            l = l-1
        end
        while x[l]-xt >= fsh
            l += 1
        end
        if xt-x[l] < fsh
            a[j] = y[l]
        else
            a[j] = freq_eval(xt, x, y, ad)
        end
    end

    nm1 = nfcns - 1   # nm1 => "nfcns minus 1"

    for j = 1 : nfcns
        dtemp = 0.0
        for k = 1 : nm1
            dtemp = muladd(a[k+1], cospi(2 * (j-1) * delf * k), dtemp)
        end
        alpha[j] = 2dtemp + a[1]
    end

    for j = 2 : nfcns
        alpha[j] *= 2.0 * delf
    end
    alpha[1] *= delf

    if !full_grid
        p[1] = muladd(2alpha[nfcns], bb, alpha[nm1])
        p[2] = 2.0*aa*alpha[nfcns]
        q[1] = alpha[nfcns-2]-alpha[nfcns]
        for j = 2 : nm1
            if j >= nm1
                aa *= 0.5
                bb *= 0.5
            end
            p[j+1] = 0.0
            for k = 1 : j
                a[k] = p[k]
                p[k] = 2.0 * bb * a[k]
            end
            p[2] = muladd(a[1], 2aa, p[2])
            for k = 1 : j-1
                p[k] += muladd(aa, a[k+1], q[k])
            end
            for k = 3 : j+1
                p[k] = muladd(aa, a[k-1], p[k])
            end

            if j != nm1
                for k = 1 : j
                    q[k] = -a[k]
                end
                q[1] += alpha[nfcns - 1 - j]
            end
        end
        for j = 1 : nfcns
            alpha[j] = p[j]
        end
    end

    if nfcns <= 3
        alpha[nfcns+1] = alpha[nfcns+2] = 0.0
    end

    #
    # CALCULATE THE IMPULSE RESPONSE.
    #
    h = zeros(Float64, nfilt)
    if !neg
        if nodd
            for j = 1 : nm1
                h[j] = 0.5 * alpha[nz-j]
            end
            h[nfcns] = alpha[1]
        else
            h[1] = 0.25 * alpha[nfcns]
            for j = 2 : nm1
                h[j] = 0.25 * (alpha[nz-j] + alpha[nfcns+2-j])
            end
            h[nfcns] = muladd(0.5, alpha[1], 0.25 * alpha[2])
        end
    else
        if nodd
            h[1] = 0.25 * alpha[nfcns]
            h[2] = 0.25 * alpha[nm1]
            for j = 1 : nm1
                h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+3-j])
            end
            h[nfcns] = muladd(0.5, alpha[1], -0.25 * alpha[3])
            h[nz] = 0.0
        else
            h[1] = 0.25 * alpha[nfcns]
            for j = 2 : nm1
                h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+2-j])
            end
            h[nfcns] = muladd(0.5, alpha[1], -0.25 * alpha[2])
        end
    end

    for j = 1 : nfcns
        k = nfilt + 1 - j
        if !neg
           h[k] = h[j]
        else
           h[k] = -h[j]
        end
    end
    if neg && nodd
        h[nz] = 0.0
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
    return remez(numtaps, band_defs; Hz=Hz, neg=neg, kwargs...)
end

