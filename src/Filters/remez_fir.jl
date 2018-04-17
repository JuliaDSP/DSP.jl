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

