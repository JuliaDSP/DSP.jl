# linear phase FIR filter design that optimizes maximum error
# in the frequency domain

"""
    remez(numtaps::Integer, bands::Array, desired::Array; 
          weight::Array=[], 
          Hz::Real=1.0, 
          filter_type::String="bandpass", 
          maxiter::Integer=25, 
          grid_density::Integer=16)

Calculate the minimax optimal filter using the Remez exchange algorithm.

Calculate the filter-coefficients for the finite impulse response
(FIR) filter whose transfer function minimizes the maximum error
between the desired gain and the realized gain in the specified
frequency bands using the Remez exchange algorithm.

Parameters
----------
numtaps : int
    The desired number of taps in the filter. The number of taps is
    the number of terms in the filter, or the filter order plus one.
bands : array_like
    A monotonic sequence containing the band edges in Hz.
    All elements must be non-negative and less than half the sampling
    frequency as given by `Hz`.
desired : array_like
    A sequence half the size of bands containing the desired gain
    in each of the specified bands.
weight : array_like, optional
    A relative weighting to give to each band region. The length of
    `weight` has to be half the length of `bands`.
Hz : scalar, optional
    The sampling frequency in Hz. Default is 1.
filter_type : {'bandpass', 'differentiator', 'hilbert'}, optional
    The type of filter:

      'bandpass' : flat response in bands. This is the default.

      'differentiator' : frequency proportional response in bands.

      'hilbert' : filter with odd symmetry, that is, type III
                  (for even order) or type IV (for odd order)
                  linear phase filters.

maxiter : int, optional
    Maximum number of iterations of the algorithm. Default is 25.
grid_density : int, optional
    Grid density. The dense grid used in `remez` is of size
    ``(numtaps + 1) * grid_density``. Default is 16.

Returns
-------
out : ndarray
    A rank-1 array containing the coefficients of the optimal
    (in a minimax sense) filter.

See Also
--------
freqz : Compute the frequency response of a digital filter.

References
----------
.. [1] J. H. McClellan and T. W. Parks, "A unified approach to the
       design of optimum FIR linear phase digital filters",
       IEEE Trans. Circuit Theory, vol. CT-20, pp. 697-701, 1973.
.. [2] J. H. McClellan, T. W. Parks and L. R. Rabiner, "A Computer
       Program for Designing Optimum FIR Linear Phase Digital
       Filters", IEEE Trans. Audio Electroacoust., vol. AU-21,
       pp. 506-525, 1973.

Examples
--------
We want to construct a filter with a passband at 0.2-0.4 Hz, and
stop bands at 0-0.1 Hz and 0.45-0.5 Hz. Note that this means that the
behavior in the frequency ranges between those bands is unspecified and
may overshoot.

>>> from scipy import signal
>>> bpass = signal.remez(72, [0, 0.1, 0.2, 0.4, 0.45, 0.5], [0, 1, 0])
>>> freq, response = signal.freqz(bpass)
>>> ampl = np.abs(response)

>>> import matplotlib.pyplot as plt
>>> fig = plt.figure()
>>> ax1 = fig.add_subplot(111)
>>> ax1.semilogy(freq/(2*np.pi), ampl, 'b-')  # freq in Hz
>>> plt.show()

"""
function remez(numtaps::Integer, bands::Array, desired::Array; 
               weight::Array=[], 
               Hz::Real=1.0, 
               filter_type::String="bandpass", 
               maxiter::Integer=25, 
               grid_density::Integer=16)

    # Convert type
    if (filter_type == "bandpass")
        tnum = 1
    elseif (filter_type == "differentiator")
        tnum = 2
    elseif (filter_type == "hilbert")
        tnum = 3
    else
        error("`filter_type` must be \"bandpass\", \"differentiator\", or \"hilbert\".")
    end
    if (length(weight)==0)
        weight = ones(desired)
    end

    h = zeros(numtaps)
    #return sigtools._remez(numtaps, bands, desired, weight, tnum, Hz,
                           #maxiter, grid_density)
    #int pre_remez(double *h2, int numtaps, int numbands, double *bands,
    #              double *response, double *weight, int type, int maxiter,
    #              int grid_density);
    bands = copy(bands)/Hz
    ccall((:pre_remez, "remez"), Integer, (Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint), h, numtaps, length(bands)/2, bands, desired, weight, tnum, maxiter, grid_density)
    return h

end



function alternating_solve(ω_grid, D_grid, extremal_indices)
    M = length(extremal_indices) - 2
    extremals = ω_grid[extremal_indices]
    
    A = zeros(Float64, M+2, M+2)
    sign_of_error = 1
    for row = 1:M+2
        for col = 1:M+1
            A[row,col] = cos((col-1)*extremals[row])
        end
        A[row,M+2] = sign_of_error
        sign_of_error *= -1
    end
    
    b = D_grid[extremal_indices]
    
    aa = A\b
    a,delta = aa[1:end-1], abs(aa[end])
    
    # inefficient calculation of DFT
    A_full = zeros(Float64, length(ω_grid), length(a))
    for row = 1:length(ω_grid)
        for col = 1:length(a)
            A_full[row,col] = cos((col-1)*ω_grid[row])
        end
    end
    H = A_full*a
    
    a, delta, H
    
end

using Formatting

function update_extremals(D, extremal_indices, band_indices_grid, H, D_grid)
    # modifies extremal_indices in-place
    for band_index = 1:length(D)   # for each band
        extremal_set = find( band_indices_grid[extremal_indices] .== band_index )
        printfmtln("band_index $band_index: {} extremals", length(extremal_set))
        for index = 2:length(extremal_set)-1  # leave the 1st and last extremal alone
            extremal_index = extremal_indices[extremal_set[index]]
            #display(extremal_index)
            #continue
            error_current = abs(H[extremal_index] - D_grid[extremal_index])
            error_left = (extremal_index>1) ? abs(H[extremal_index-1] - D_grid[extremal_index-1]) : error_current
            error_right = (extremal_index<length(H)) ? abs(H[extremal_index+1] - D_grid[extremal_index+1]) : error_current
            if error_left > error_current
                increment = -1
                error_next = error_left
            elseif error_right > error_current
                increment = 1
                error_next = error_right
            else
                increment = 0
                error_next = error_current
            end
            extremal_index += increment
            while error_next > error_current
                error_current = error_next
                extremal_index += increment
                index_next = extremal_index + increment
                if (index_next < 1) || (index_next > length(D_grid))
                    continue
                end
                error_next = abs(H[extremal_index+increment] - D_grid[extremal_index+increment])
            end
            #display([extremal_indices[extremal_set[index]] extremal_index])
            extremal_indices[extremal_set[index]] = extremal_index
        end
    end
    #display(extremal_indices)
end


function remez_jl(numtaps::Integer, bands::Array, desired::Array; 
                  weight::Array=[], 
                  Hz::Real=1.0, 
                  maxiter::Integer=25, 
                  grid_density::Integer=16)

    if (length(weight)==0)
        weight = ones(desired)
    end

    bands = copy(bands)/Hz

    
    L = numtaps
    M = Int(floor(L / 2))
    grid_spacing = π / (grid_density * M)
    ω_pairs = 2π * bands
    D = desired
    
    ω_grid = zeros(Float64, 0)
    D_grid = zeros(Float64, 0)
    band_indices_grid = zeros(Int64, 0)
    band_edge_indices = zeros(Int64, 0)
    for band_index = 1:length(D)
        push!(band_edge_indices, length(ω_grid)+1)
        pair_index = (band_index-1)*2 + 1
        a, b = ω_pairs[pair_index:pair_index+1]
        D_band = D[band_index]
        ω_grid_band = collect(a:grid_spacing:b)
        if (ω_grid_band[end] < b)
            push!(ω_grid_band, b)
        end
        append!(ω_grid, ω_grid_band)
        append!(D_grid, D_band*(1+0*ω_grid_band))
        append!(band_indices_grid, band_index*ones(Int64, length(ω_grid_band)))
        push!(band_edge_indices, length(ω_grid))
    end

    # Initial extremal frequencies
    extremal_indices = Array{Int,1}(round.(collect(linspace(1,length(ω_grid),M+2))))  # evenly spaced
    extremal_bands = band_indices_grid[extremal_indices]
    # now add band edges to set
    for band_index = 1:length(D)
        extremal_set = find( band_indices_grid[extremal_indices] .== band_index )
        extremal_indices_in_band = extremal_indices[extremal_set]
        pair_index = 2*(band_index-1) + 1
        extremal_indices[extremal_set[1]] = band_edge_indices[pair_index]
        extremal_indices[extremal_set[end]] = band_edge_indices[pair_index+1]
    end
    #extremals = ω_grid[extremal_indices]

    a, delta, H = alternating_solve(ω_grid, D_grid, extremal_indices)
    println("delta $delta")

    iteration = 0
    while iteration < maxiter
        println("iteration $iteration")
        update_extremals(D, extremal_indices, band_indices_grid, H, D_grid)
        a, delta, H = alternating_solve(ω_grid, D_grid, extremal_indices)
        println("delta $delta")
        iteration += 1
    end
    
    a
end






function lagrange_interp(k::Integer, n::Integer, m::Integer, x::AbstractVector)
    retval = 1.0;
    q = x[k];
    for l = 1 : m
        for j = l : m : n
            if (j != k)
                retval *= 2.0 * (q - x[j])
            end
        end
    end
    1.0 / retval
end


"""
eff(freq::Float64, fx::AbstractVector, lband::Integer, jtype::Integer)

/*
 *-----------------------------------------------------------------------
 * FUNCTION: eff
 *  FUNCTION TO CALCULATE THE DESIRED MAGNITUDE RESPONSE
 *  AS A FUNCTION OF FREQUENCY.
 *  AN ARBITRARY FUNCTION OF FREQUENCY CAN BE
 *  APPROXIMATED IF THE USER REPLACES THIS FUNCTION
 *  WITH THE APPROPRIATE CODE TO EVALUATE THE IDEAL
 *  MAGNITUDE.  NOTE THAT THE PARAMETER FREQ IS THE
 *  VALUE OF NORMALIZED FREQUENCY NEEDED FOR EVALUATION.
 *-----------------------------------------------------------------------
 */
"""
function eff(freq::Float64, fx::AbstractVector, lband::Integer, jtype::Integer)
    if (jtype != 2)
        return fx[lband]
    else 
        return fx[lband] * freq
    end
end

"""
wate(freq::Float64, fx::AbstractVector, wtx::AbstractVector, lband::Integer, jtype::Integer)
/*
 *-----------------------------------------------------------------------
 * FUNCTION: wate
 *  FUNCTION TO CALCULATE THE WEIGHT FUNCTION AS A FUNCTION
 *  OF FREQUENCY.  SIMILAR TO THE FUNCTION eff, THIS FUNCTION CAN
 *  BE REPLACED BY A USER-WRITTEN ROUTINE TO CALCULATE ANY
 *  DESIRED WEIGHTING FUNCTION.
 *-----------------------------------------------------------------------
 */
"""
function wate(freq::Float64, fx::AbstractVector, wtx::AbstractVector, lband::Integer, jtype::Integer)
    if (jtype != 2)
        return wtx[lband]
    end
    if (fx[lband] >= 0.0001)
        return wtx[lband] / freq
    end
    wtx[lband]
end


"""
build_grid(numtaps, bands, desired, weight, grid_density)
return "grid" and "des" and "wt" arrays
"""
function build_grid(numtaps, bands, desired, weight, grid_density)
    # translate from scipy remez argument names
    L = numtaps
    M = Int(floor(L / 2))
    grid_spacing = 0.5 / (grid_density * M)  # "delf" or "delta f"

    lgrid = grid_density
    dimsize = Int( ceil(numtaps/2.0 + 2) )
    wrksize = grid_density * dimsize

    grid = zeros(Float64, wrksize)  # the array of frequencies, between 0 and 0.5
    des = zeros(Float64, wrksize)   # the desired function on the grid
    wt = zeros(Float64, wrksize)    # array of weights

    fx = desired
    wtx = weight
    
    jtype = 1     # FIXME should be an input
    nbands = length(desired)
    edge = bands
    nfilt = numtaps

    neg = 1
    if (jtype == 1)
        neg = 0
    end
    nodd = nfilt % 2
    nfcns = nfilt ÷ 2
    if (nodd == 1 && neg == 0) 
        nfcns = nfcns + 1
    end

    """/*
     * SET UP THE DENSE GRID. THE NUMBER OF POINTS IN THE GRID
     * IS (FILTER LENGTH + 1)*GRID DENSITY/2
     */"""
    grid[1] = edge[1]
    delf = lgrid * nfcns
    delf = 0.5 / delf
    if (neg != 0)
        if (edge[1] < delf)
            grid[1] = delf
        end
    end
    j = 1
    l = 1
    lband = 1

    """/*
     * CALCULATE THE DESIRED MAGNITUDE RESPONSE AND THE WEIGHT
     * FUNCTION ON THE GRID
     */"""
    while true
        fup = edge[l + 1]
        while true
            temp = grid[j]
            des[j] = eff(temp,fx,lband,jtype)
            wt[j] = wate(temp,fx,wtx,lband,jtype)
            j += 1
            if (j > wrksize)
                # too many points, or too dense grid
                return -1
            end
            grid[j] = temp + delf
            if (grid[j] > fup)
                break
            end
        end
        grid[j-1] = fup
        des[j-1] = eff(fup,fx,lband,jtype)
        wt[j-1] = wate(fup,fx,wtx,lband,jtype)
        lband += 1
        l += 2
        if (lband > nbands) 
            break;
        end
        grid[j] = edge[l]
    end

    ngrid = j - 1
    if (neg == nodd)
        if (grid[ngrid] > (0.5-delf))
            ngrid -= 1
        end
    end
    
    grid[1:ngrid], des[1:ngrid], wt[1:ngrid]

end

"""
function freq_eval(k::Integer, n::Integer, grid::AbstractVector, 
                   x::AbstractVector, y::AbstractVector, ad::AbstractVector)
-----------------------------------------------------------------------
 FUNCTION: freq_eval (gee)
  FUNCTION TO EVALUATE THE FREQUENCY RESPONSE USING THE
  LAGRANGE INTERPOLATION FORMULA IN THE BARYCENTRIC FORM
-----------------------------------------------------------------------
"""
function freq_eval(k::Integer, n::Integer, grid::AbstractVector, x::AbstractVector, y::AbstractVector, ad::AbstractVector)
    d = 0.0;
    p = 0.0;
    xf = cos((2π)*grid[k])

    for j = 1 : n
        c = ad[j] / (xf - x[j])
        d += c
        p += c * y[j]
    end

    p/d
end

function initialize_y(dev::Float64, nz::Integer, iext::AbstractArray, des::AbstractArray, wt::AbstractArray, y::AbstractArray)
    nu = 1
    if ( (dev) > 0.0 )
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


mutable struct ExtremalSet
    j::Int  # index of the current extremal being updated
    nz::Int # number of cosines in the approximation (including the constant term).
            # nz = nfcns + 1 where nfcns = nfilt / 2, and 
            # nfilt is the filter length or number of taps. 
            # For example, for a length 15 filter, nfcns = 7 and nz = 8. 
    jchgne::Int  # number of extremal indices that changed this iteration
    k1::Int
    knz::Int
    klow::Int
    nu::Int
    nut::Int
    l::Int
    luck::Int
    iext::Array{Int64}  # the list of extremal indices
end

struct Goto
    L200
    L210
    L215
    L220
    L225
    L230
    L235
    L240
    L250
    L260
    L300
    L310
    L315
    L320
    L325
    L330
    L340
    L350
    L370
end





function remez_jl2(numtaps::Integer, bands::Array, desired::Array; 
                   weight::Array=[], 
                   Hz::Real=1.0, 
                   maxiter::Integer=25, 
                   grid_density::Integer=16)

    if (length(weight)==0)
        weight = ones(desired)
    end
    bands = copy(bands)/Hz

    bands = vec(bands)
    desired = vec(desired)
    weight = vec(weight)
    
    grid, des, wt = build_grid(numtaps, bands, desired, weight, grid_density);

    ngrid = length(grid)
    
    #for j = 1 : ngrid
    #    printfmtln("  j={}: grid[j]={}", j, grid[j]);
    #end
    
    jtype = 1   # should be input, and passed down to remez
    neg = 1;
    if (jtype == 1)
      neg = 0
    end
    nodd = numtaps % 2
    nfcns = numtaps ÷ 2  # integer divide
    if (nodd == 1 && neg == 0)
        nfcns = nfcns + 1
    end
    temp = Float64(ngrid-1) / nfcns
    dimsize = Int64(ceil(numtaps/2.0 + 2))
    printfmtln("  dimsize={}   nfcns={}", dimsize, nfcns);
    iext = vec(zeros(Int64, dimsize))   # indices of extremals
    for j = 1 : nfcns
        iext[j] = Int64(floor((j-1)*temp)) + 1
    end
    iext[nfcns+1] = ngrid;

    devl = -1.0;  # "dev" deviation on last iteration
    dev = 0.0
    nz  = nfcns+1;
    nzz = nfcns+2;
    niter = 0;
    
    jet = ((nfcns-1) ÷ 15) + 1

    while true
    
        #
        # Start next iteration
        #        
        iext[nzz] = ngrid + 1;
        niter += 1;
        if (niter > maxiter) 
            break
        end
        
        x = cos.(grid[iext[1:nz]]*2*π)
        for j = 1 : nz
            printfmtln("  j={}: iext[j]={} grid[iext[j]]={} x[j]={}", j, iext[j], grid[iext[j]], x[j]);
        end
        
        nz = nfcns + 1
        ad = zeros(Float64, nz)
        for j = 1 : nz
            ad[j] = lagrange_interp(j, nz, jet, x)
        end

        dnum = 0.0
        dden = 0.0
        k = 1
        for j = 1 : nz
            l = iext[j]
            dnum += ad[j] * des[l]
            dden += k * ad[j] / wt[l]
            k = -k
        end
        dev = dnum / dden
        printfmtln("DEVIATION = {}", dev);

        y = 0*x
        nu, dev = initialize_y(dev, nz, iext, des, wt, y)
        
        if ( dev <= devl )
            # finished 
            return -1  # error - deviation should always increase
        end
        
        devl = dev;
        jchnge = 0;
        k1 = iext[1];
        knz = iext[nz];
        klow = 0;
        nut = -nu;
        j = 1;

        #
        # SEARCH FOR THE EXTREMAL FREQUENCIES OF THE BEST APPROXIMATION
        #

        if (jchnge <= 0)  # we are done if none of the extremal indices changed
            break
        end
    end  # while
    return dev
end



