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




using Formatting




"""/*
 *-----------------------------------------------------------------------
 * FUNCTION: lagrange_interp (d)
 *  FUNCTION TO CALCULATE THE LAGRANGE INTERPOLATION
 *  COEFFICIENTS FOR USE IN THE FUNCTION gee.
 *-----------------------------------------------------------------------
 */"""
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
function build_grid(numtaps, bands, desired, weight, grid_density, jtype)
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
    
    #jtype = 1
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
function remez_jl2(numtaps::Integer, bands::Array, desired::Array; 
                   weight::Array=[], 
                   Hz::Real=1.0, 
                   filter_type::String="bandpass",
                   maxiter::Integer=25, 
                   grid_density::Integer=16)
    # Convert type
    if (filter_type == "bandpass")
        jtype = 1
    elseif (filter_type == "differentiator")
        jtype = 2
    elseif (filter_type == "hilbert")
        jtype = 3
    else
        error("`filter_type` must be \"bandpass\", \"differentiator\", or \"hilbert\".")
    end

    if (length(weight)==0)
        weight = ones(desired)
    end
    bands = copy(bands)/Hz

    bands = vec(bands)   # in C, known as "edge"
    desired = vec(desired)
    weight = vec(weight)
    
    grid, des, wt = build_grid(numtaps, bands, desired, weight, grid_density, jtype);

    nfilt = numtaps
    ngrid = length(grid)
    
    #for j = 1 : ngrid
    #    printfmtln("  j={}: grid[j]={}", j, grid[j]);
    #end
    
    # "jtype" is the type of filter, with the following meaning.
      #define BANDPASS       1
      #define DIFFERENTIATOR 2
      #define HILBERT        3
    # I think the "j" is because it is an int, and FORTRAN assigned
    # types based on starting letter of variable name LOL.
    # jtype input:
    #    Type I and II symmetric linear phase: neg==0   (jtype==1)
    #    Type III and IV negative symmetric linear phase: neg==1   (jtype==2 or 3)
    neg = 1     # "neg" means negative symmetry.
    if (jtype == 1)
      neg = 0
    end
    nodd = numtaps % 2   # nodd == 1: length is odd; nodd == 0: length is even
    nfcns = numtaps ÷ 2  # integer divide
    if (nodd == 1 && neg == 0)
        nfcns = nfcns + 1
    end
    temp = Float64(ngrid-1) / nfcns
    dimsize = Int64(ceil(numtaps/2.0 + 2))
    #printfmtln("  dimsize={}   nfcns={}", dimsize, nfcns);

    """/*
     * SET UP A NEW APPROXIMATION PROBLEM WHICH IS EQUIVALENT
     * TO THE ORIGINAL PROBLEM
     */"""
    if (neg <= 0)
        if (nodd != 1)
            for j = 1 : ngrid
                change = cos(π*grid[j]);
                des[j] = des[j] / change;
                wt[j]  = wt[j] * change;
            end
        end
    else
        if (nodd != 1)
            for j = 1 : ngrid
                change = sin(π*grid[j]);
                des[j] = des[j] / change;
                wt[j]  = wt[j]  * change;
            end
        else
            for j = 1 : ngrid
                change = sin(2π * grid[j]);
                des[j] = des[j] / change;
                wt[j]  = wt[j]  * change;
            end
        end
    end
        
    iext = vec(zeros(Int64, dimsize))   # indices of extremals
    x = vec(zeros(Float64, dimsize))
    y = vec(zeros(Float64, dimsize))

    for j = 1 : nfcns
        iext[j] = Int64(floor((j-1)*temp)) + 1
    end
    iext[nfcns+1] = ngrid;

    dev = 0.0     # deviation from the desired function, 
                  # that is, the amount of "ripple" on the extremal set
    devl = -1.0;  # deviation on last iteration
    nz  = nfcns+1;
    nzz = nfcns+2;
    niter = 0;
    ad = zeros(Float64, nz)
    
    jet = ((nfcns-1) ÷ 15) + 1

    while true
    
        #
        # Start next iteration
        #        
      @label L100
        iext[nzz] = ngrid + 1;
        niter += 1;
        if (niter > maxiter) 
            break
        end
        
        x[1:nz] = cos.(grid[iext[1:nz]]*2*π)
        #for j = 1 : nz
        #    printfmtln("  j={}: iext[j]={} grid[iext[j]]={} x[j]={}", j, iext[j], grid[iext[j]], x[j]);
        #end
        
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
        # printfmtln("DEVIATION = {}", dev);

        y = 0*x
        nu, dev = initialize_y(dev, nz, iext, des, wt, y)
        
        if ( dev <= devl )
            # finished 
            return -1  # error - deviation should always increase
        end        
        devl = dev;

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
        jchnge = 0;
        k1 = iext[1];
        knz = iext[nz];
        klow = 0;
        nut = -nu;
        j = 1;

      @label L200
        if (j == nzz) ynz = comp; end
        if (j >= nzz) @goto L300; end
        kup = iext[j+1];
        l = iext[j]+1;
        nut = -nut;
        if (j == 2) y1 = comp; end
        comp = dev;
        if (l >= kup) @goto L220; end
        err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
        if ((nut*err-comp) <= 0.0) @goto L220; end
        comp = nut * err;
      @label L210
        l += 1; if (l >= kup) @goto L215; end
        err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
        if ((nut*err-comp) <= 0.0) @goto L215; end
        comp = nut * err;
        @goto L210;

      @label L215
        iext[j] = l - 1; j += 1
        klow = l - 1;
        jchnge += 1;
        @goto L200;

      @label L220
        l -= 1;
      @label L225
        l -= 1; if (l <= klow) @goto L250; end
        err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
        if ((nut*err-comp) > 0.0) @goto L230; end
        if (jchnge <= 0) @goto L225; end
        @goto L260;

      @label L230
        comp = nut * err;
      @label L235
        l -= 1; if (l <= klow) @goto L240; end
        err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
        if ((nut*err-comp) <= 0.0) @goto L240; end
        comp = nut * err;
        @goto L235;
      @label L240
        klow = iext[j];
        iext[j] = l+1;
        j += 1;
        jchnge += 1;
        @goto L200;

      @label L250
        l = iext[j]+1;
        if (jchnge > 0) @goto L215; end

      @label L255
        l += 1; if (l >= kup) @goto L260; end
        err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
        if ((nut*err-comp) <= 0.0) @goto L255; end
        comp = nut * err;

        @goto L210;
      @label L260
        klow = iext[j]; j += 1
        @goto L200;

      @label L300
        if (j > nzz) @goto L320; end
        if (k1 > iext[1] ) k1 = iext[1]; end
        if (knz < iext[nz]) knz = iext[nz]; end
        nut1 = nut;
        nut = -nu;
        l = 0;
        kup = k1;
        comp = ynz*(1.00001);
        luck = 1;
      @label L310
        l += 1; if (l >= kup) @goto L315; end
        err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
        if ((nut*err-comp) <= 0.0) @goto L310; end
        comp =  nut * err;
        j = nzz;
        @goto L210;

      @label L315
        luck = 6;
        @goto L325;

      @label L320
        if (luck > 9) @goto L350; end
        if (comp > y1) y1 = comp; end
        k1 = iext[nzz];
      @label L325
        l = ngrid+1;
        klow = knz;
        nut = -nut1;
        comp = y1*(1.00001);
      @label L330
        l -= 1; if (l <= klow) @goto L340; end
        err = (freq_eval(l,nz,grid,x,y,ad)-des[l]) * wt[l];
        if ((nut*err-comp) <= 0.0) @goto L330; end
        j = nzz;
        comp =  nut * err;
        luck = luck + 10;
        @goto L235;
      @label L340
        if (luck == 6) @goto L370; end
        for j = 1 : nfcns
            iext[nzz-j] = iext[nz-j];
        end
        iext[1] = k1;
        @goto L100;
      @label L350
        kn = iext[nzz];
        for j = 1 : nfcns
            iext[j] = iext[j+1];
        end
        iext[nz] = kn;

        @goto L100;
      @label L370
        ;
        
        if (jchnge <= 0)  # we are done if none of the extremal indices changed
            break
        end
    end  # while
    
    # 
    #    CALCULATION OF THE COEFFICIENTS OF THE BEST APPROXIMATION
    #    USING THE INVERSE DISCRETE FOURIER TRANSFORM
    #

    a = vec(zeros(Float64, dimsize))   # frequency response on evenly spaced grid
    p = vec(zeros(Float64, dimsize))
    q = vec(zeros(Float64, dimsize)) 
    alpha = vec(zeros(Float64, dimsize))   # return vector
    
    PI = 3.14159265358979323846
    TWOPI = (PI+PI)
    nm1 = nfcns - 1;   # nm1 => "nfcns minus 1"
    fsh = 1.0e-06;
    gtemp = grid[1];   # grid[1] is temporarily used in freq_eval in this loop
    x[nzz] = -2.0;
    cn  = 2*nfcns - 1;
    delf = 1.0/cn;
    l = 1;
    kkk = 0;

    if (bands[1] == 0.0 && bands[end] == 0.5) kkk = 1; end

    if (nfcns <= 3) kkk = 1; end
    if (kkk !=     1)
        dtemp = cos(TWOPI*grid[1]);
        dnum  = cos(TWOPI*grid[ngrid]);
        aa    = 2.0/(dtemp-dnum);
        bb    = -(dtemp+dnum)/(dtemp-dnum);
    end

    # Fill in "a" array with the frequency response on an evenly
    # spaced set of frequencies. Care is taken to use "y[l]" -
    # an already computed response on one of the extremals "l" -
    # if the extremal is equal to the frequency ft. If no y[l]
    # matches, a[j] is computed using freq_eval.
    for j = 1 : nfcns 
        ft = (j - 1) * delf;
        xt = cos(TWOPI*ft);
        if (kkk != 1)
            xt = (xt-bb)/aa;
            ft = acos(xt)/TWOPI;
        end
      @label L410
        xe = x[l];
        if (xt > xe) @goto L420; end
        if ((xe-xt) < fsh) @goto L415; end
        l += 1;
        @goto L410;
      @label L415
        a[j] = y[l];
        @goto L425;
      @label L420
        if ((xt-xe) < fsh) @goto L415; end
        grid[1] = ft;
        a[j] = freq_eval(1,nz,grid,x,y,ad);
      @label L425
        if (l > 1) l = l-1; end
    end
    grid[1] = gtemp;  # restore grid[1]
    #for j = 1 : nfcns 
    #    printfmtln("  j={}: a[j]={}", j, a[j]);
    #end

    dden = TWOPI / cn;
    for j = 1 : nfcns
        dtemp = 0.0;
        dnum = (j-1) * dden;
        if (nm1 >= 1)
            for k = 1 : nm1
                dtemp += a[k+1] * cos(dnum*k);
            end
        end
        alpha[j] = 2.0 * dtemp + a[1];
    end

    for j = 2 : nfcns
        alpha[j] *= 2.0 / cn;
    end
    alpha[1] /= cn;

    if (kkk != 1)
        p[1] = 2.0*alpha[nfcns]*bb+alpha[nm1];
        p[2] = 2.0*aa*alpha[nfcns];
        q[1] = alpha[nfcns-2]-alpha[nfcns];
        for j = 2 : nm1
            if (j >= nm1)
                aa *= 0.5;
                bb *= 0.5;
            end
            p[j+1] = 0.0;
            for k = 1 : j
                a[k] = p[k];
                p[k] = 2.0 * bb * a[k];
            end
            p[2] += a[1] * 2.0 *aa;
            jm1 = j - 1;
            for k = 1 : jm1
                p[k] += q[k] + aa * a[k+1];
            end
            jp1 = j + 1;
            for k = 3 : jp1
                p[k] += aa * a[k-1];
            end

            if (j != nm1)
                for k = 1 : j
                    q[k] = -a[k];
                end
                q[1] += alpha[nfcns - 1 - j];
            end
        end
        for j = 1 : nfcns 
            alpha[j] = p[j];
        end
    end

    if (nfcns <= 3)
        alpha[nfcns+1] = alpha[nfcns+2] = 0.0;
    end

    #for j = 1 : nfcns 
    #    printfmtln("  j={}: alpha[j]={}", j, alpha[j]);
    #end
    
    #
    # CALCULATE THE IMPULSE RESPONSE.
    #
    h = zeros(Float64, nfilt)
    if (neg <= 0)
        if (nodd != 0)
            for j = 1 : nm1
                h[j] = 0.5 * alpha[nz-j];
            end
            h[nfcns] = alpha[1];
        else
            h[1] = 0.25 * alpha[nfcns];
            for j = 2 : nm1
                h[j] = 0.25 * (alpha[nz-j] + alpha[nfcns+2-j]);
            end
            h[nfcns] = 0.5*alpha[1] + 0.25*alpha[2];
        end
    else
        if (nodd != 0)
            h[1] = 0.25 * alpha[nfcns];
            h[2] = 0.25 * alpha[nm1];
            for j = 1 : nm1
                h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+3-j]);
            end
            h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[3];
            h[nz] = 0.0;
        else
            h[1] = 0.25 * alpha[nfcns];
            for j = 2 : nm1
                h[j] = 0.25 * (alpha[nz-j] - alpha[nfcns+2-j]);
            end
            h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[2];
        end
    end

    for j = 1 : nfcns
        k = nfilt + 1 - j;
        if (neg == 0)
           h[k] = h[j];
        else
           h[k] = -h[j];
        end
      end
    if (neg == 1 && nodd == 1) h[nz] = 0.0; end
    
    return h
end



