DSP.jl provides a number of common DSP routines in Julia.  So far, the following functions are
implemented:

digital filtering:  
* filt

correlation and convolution:  
* conv  
* conv2  
* deconv 
* xcorr 

FFTs provided by FFTW interface:
* bfft
* bfftn
* brfft
* brfftn
* fft
* fft2
* fft3
* fftn
* ifft
* ifft2
* ifft3
* ifftn
* irfft
* irfftn
* rfft
* rfftn

FFT utilities:
* fftshift	
* ifftshift

periodogram estimation:
* periodogram
* welch_pgram
* bartlett_pgram

window functions:
* rect
* hanning
* hamming
* tukey
* cosine
* lanczos
* triang
* bartlett
* gaussian
* bartlett_hann
* blackman
* kaiser

common DSP mathematics:
* sinc

auxiliary functions:
* arraysplit	
