
module Periodograms

export Periodogramt, power, freq, spectrogram, time, pf, stf, periodogram, welch_pgram, bartlett_pgram

typealias DSPNumber FFTW.fftwNumber
typealias DSPReal FFTW.fftwReal
typealias DSPComplex FFTW.fftwComplex

type Periodogramt{T<:DSPNumber}   # Periodogram is the module name, so there is a conflict
    x::Array{T,1}
    n::Integer                      # 0 < n <= length(s)
    noverlap::Integer               # 0 <= noverlap < n
    fs::Real                        # positive
    window::Union(Function,Bool)    # false || Function
    twosided::Bool                  # false for 1-sided
    taper::Real                     # 0 <= taper <=0.5, taper edges for periodograms or when using pad, int(n*taper) elements are tapered on each side using split cosine bell
    #pad::Integer                   # 0 <= pad < n, zero padding of windows
                                    # takes a section from x of length n-pad, 
                                    # and pads with zero to give padded section of length n
end
function Periodogramt{T<:DSPNumber}(x::Array{T,1}; n = int(length(x)/8), noverlap = int(length(x)/16), fs = 1, window = DSP.Windows.hamming, twosided = false, taper = 0)
    eltype(x)<:DSPComplex && (twosided=true)  # make it 2-sided for complex data
    po = Periodogramt(x,n,noverlap,fs,window,twosided,taper)
    errorcheck(po)
    return po
end
function Periodogramt{T<:DSPNumber}(x::Array{T,2}; args...)
    if size(x,1)==1 || size(x,2)==1
        return Periodogramt(vec(x);args...)  # convert to a vector
    else
        error("not supported for 2D arrays")
    end
end
function errorcheck{T<:DSPNumber}(P::Periodogramt{T})
    !(0 < P.n <= length(P.x)) && error("use 0 < n <= length(x)")
    !(0 <= P.noverlap < P.n) && error("use 0 <= noverlap < n")
    P.fs <= 0 && error("fs is not positive")
    P.window == true && error("window = true not supported")
    eltype(P.x)<:DSPComplex && P.twosided==false && error("twosided=true for complex signal not supported")
end

# ======= the old DSP functions =======

# sided=2 default for now since that is the current implementation in DSP
# could be changed to call directly power(Periodogramt...)
function periodogram{T<:DSPNumber}(s::Array{T,1})
    Periodogramt(s,n=length(s),noverlap=0,twosided=true,window=false)
end
function periodogram{T<:DSPNumber}(s::Array{T,1}, window)
    Periodogramt(s,n=length(s),noverlap=0,twosided=true,window=window)
end
function welch_pgram{T<:DSPNumber}(s::Array{T,1}, n, m )
    Periodogramt(s,n=n,noverlap=m,twosided=true,window=false)
end
function welch_pgram{T<:DSPNumber}(s::Array{T,1}, n, m, window )
    Periodogramt(s,n=n,noverlap=m,window=window,twosided=true)
end
function bartlett_pgram{T<:DSPNumber}(s::Array{T,1}, n)
    welch_pgram(s, n, 0)
end
function bartlett_pgram{T<:DSPNumber}(s::Array{T,1}, n, window)
    welch_pgram(s, n, 0, window)
end
function spectrogram{T<:DSPNumber}(s::Array{T,1}; n=int(length(s)/8), m=int(n/2), r=1, w=(n)->ones(n))
    Periodogramt(s,n=n,noverlap=m,fs=r,window=w,twosided=true)
end

# ======= methods to apply to type object =======

# spectrum vector
function power{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    return power_w(P.x,P.n,P.noverlap,P.fs,P.window,P.twosided)
end
# frequency vector
function freq{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    n ,fs = P.n, P.fs
    if !P.twosided
        nf = div(n,2)+1 # the highest frequency + 1
    else
        nf = n
    end
    f = Array(typeof(1.0),nf)
    fsn = fs/n
    for i = 1:length(f)
        if i <= div(n,2)+1
            @inbounds f[i] = (i-1)*fsn
        else
            @inbounds f[i] = (-N+i-1)*fsn
        end
    end
    return f
end
# spectrogram matrix
function spectrogram{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    return sp_gram(P.x,P.n,P.noverlap,P.fs,P.window,P.twosided)
end
# time vector (for spectrogram)
function time{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    index = arrayspliti(P.x, P.n, P.noverlap)
    return ( (0:length(index)-1)*(P.n-P.noverlap) + P.n/2) / P.fs
end
# p, f = pf(P) for power and frequency in one command
function pf{T<:DSPNumber}(P::Periodogramt{T})
    return power(P), freq(P)
end
# s, t, f = stf(P) for spectrogram, time and frequency in one command
function stf{T<:DSPNumber}(P::Periodogramt{T})
    return spectrogram(P), time(P), freq(P)
end

# ======= implementation of power and freq methods =======
# ======= not to be exported?

# welch type spectrum (or periodogram type depending on parameters)
function power_w{T<:DSPNumber}(s::Array{T,1}, n::Integer, m::Integer, r::Real, window::Union(Function,Bool), twosided::Bool)
    tsReal = typeof(real(s[1]))
    tsComplex = typeof(complex(s[1]))
    if isa(window,Function)
        w = window(n)
    end
    index = arrayspliti(s, n, m)
    welch = false
    length(index)>1 && (welch = true)
    if twosided
        if welch
            p = zeros(tsReal,n)
        end
        split = Array(tsComplex,n)
        for i in index
            if isa(window,Function)
                copy!(split,1,s,i,n)
                broadcast!(*, split, split, w)
            else
                copy!(split,1,s,i,n)
            end
            periodogramtl!(split)
            if welch
                p += real(split)
            else
                p = real(split)
            end
        end
    else # 1-sided
        np = div(n,2)+1
        if welch
            p = zeros(tsReal,np)
        end
        split = Array(tsReal,n)
        splitrf = Array(tsComplex, np)
        for i in index
            if isa(window,Function)
                copy!(split,1,s,i,n)
                broadcast!(*, split, split, w)
            else
                copy!(split,1,s,i,n)
            end
            periodogramtl1s!(splitrf, split)
            if welch
                p += real(splitrf)
            else
                p = real(splitrf)
            end
        end
    end
    if isa(window,Function)
        p /= (r*length(index)*sum(w.^2))
    else
        p /= (r*length(index)*n)
    end
    return p
end

# 2 sided
function periodogramtl!{T<:DSPComplex}(s::Array{T,1})
    fft!(s)
    s[:] = abs2(s)
    return nothing
end
# 1 sided
function periodogramtl1s!{T<:DSPComplex,T2<:DSPReal}(sr::Array{T,1},s::Array{T2,1})
    sr[:] = rfft(s)
    sr[:] = abs2(sr)*2
    sr[1] /= 2
    iseven(length(s)) && (sr[end]/=2)
    return nothing
end
# index for array split
function arrayspliti{T<:DSPNumber}(s::Array{T,1}, n::Integer, m::Integer)
    l = n - m
    k = ifloor(length(s)/l - n/l + 1)
    [a*l + 1 for a=0:(k-1)]
end

# spectrogram
function sp_gram{T<:DSPNumber}(s::Array{T,1}, n::Integer, m::Integer, r::Real, window::Union(Function,Bool), twosided::Bool)
    tsReal = typeof(real(s[1]))
    tsComplex = typeof(complex(s[1]))
    if isa(window,Function)
        w = window(n)
    end
    index = arrayspliti(s, n, m)
    if twosided
        p = Array(tsReal,n,length(index))
        split = Array(tsComplex,n)
        j = 1
        for i in index
            if isa(window,Function)
                copy!(split,1,s,i,n)
                broadcast!(*, split, split, w)
            else
                copy!(split,1,s,i,n)
            end
            periodogramtl!(split)
            p[:,j] = real(split)
            j += 1
        end
    else # 1-sided
        np = div(n,2)+1
        p = zeros(tsReal,np,length(index))
        split = Array(tsReal,n)
        splitrf = Array(tsComplex, np)
        j = 1
        for i in index
            if isa(window,Function)
                copy!(split,1,s,i,n)
                broadcast!(*, split, split, w)
            else
                copy!(split,1,s,i,n)
            end
            periodogramtl1s!(splitrf, split)
            p[:,j] = real(splitrf)
            j += 1
        end
    end
    if isa(window,Function)
        p /= (r*sum(w.^2))
    else
        p /= (r*n)
    end
    return p
end

# split cosine bell taper window
function splitcosinebell(m::Integer)
    return [0.5*(1-cos(pi*(t+0.5)/m)) for t in 0:m-1]
end
# apply taper of length mt on both sides of s[1:nt], s[nt+1:end] is unchanged
function applytaper!{T<:DSPNumber}(s::Array{T,1}, nt::Integer, mt::Integer)
    int(mt*2)>nt && error("mt too big")
    length(s)<nt && error("nt too big")
    scb = splitcosinebell(mt)
    for i = 1:mt
        @inbounds s[i] *= scb[i]
        @inbounds s[nt-i+1] *= scb[i]
    end
end


# =======================================================


# testing
using Base.Test

# timing of methods 
using DSP
println("\n welch periodograms : \n")
N=1024*32
x0=rand(N)
nn=div(N,8)
mm=div(N,16)
P=Periodogramt(x0)  #default everything
Pp = power(P)
# initialize variables before timing
P=Periodogramt(x0,n=nn,noverlap=mm,twosided=true,window=n->ones(n))
Pp = power(P)
P.twosided=false
Pp1 = power(P)
P.twosided=true
wp = DSP.Periodogram.welch_pgram(x0,nn,mm)

tn=10

println("  P 2 sided")
@time begin for i=1:tn
Pp = power(P)
end
end
println(typeof(Pp)," ",length(Pp))

println("  old welch")
@time begin for i=1:tn
wp = DSP.Periodogram.welch_pgram(x0,nn,mm)
end
end
println(typeof(wp)," ",length(wp))

println("  P 1 sided")
P.twosided=false
@time begin for i=1:tn
Pp1 = power(P)
end
end
println(typeof(Pp1)," ",length(Pp1))

@test_approx_eq Pp wp

println("\n periodograms : \n")
nn=N
mm=0
P=Periodogramt(x0,n=nn,noverlap=mm,twosided=true,window=false)
Pp = power(P)
P.twosided=false
Pp1 = power(P)
P.twosided=true
wp = DSP.Periodogram.periodogram(x0)

println("  P 2 sided")
@time begin for i=1:tn
Pp = power(P)
end
end
println(typeof(Pp)," ",length(Pp))

println("  old periodogram")
@time begin for i=1:tn
wp = DSP.Periodogram.periodogram(x0)
end
end
println(typeof(wp)," ",length(wp))

println("  P 1 sided")
P.twosided=false
@time begin for i=1:tn
Pp1 = power(P)
end
end
println(typeof(Pp1)," ",length(Pp1))

@test_approx_eq Pp wp


println("\n spectrograms : \n")

nn=div(N,8)
mm=div(N,16)
# initialize variables before timing
P=Periodogramt(x0,n=nn,noverlap=mm,twosided=true,window=n->ones(n))
Ps = spectrogram(P)
sp,t,f = DSP.Periodogram.spectrogram(x0,n=nn,m=mm)
P.twosided=false
Ps1 = spectrogram(P)

println("  P spectrogram")
P.twosided=true
@time begin for i=1:tn
Ps = spectrogram(P)
end
end
println(typeof(Ps)," ",length(Ps))

println("  old spectrogram")
@time begin for i=1:tn
sp,t,f = DSP.Periodogram.spectrogram(x0,n=nn,m=mm)
end
end
println(typeof(sp)," ",length(sp))

println("  P spectrogram 1 sided")
P.twosided=false
@time begin for i=1:tn
Ps1 = spectrogram(P)
end
end
println(typeof(Ps1)," ",length(Ps1))

@test_approx_eq Ps sp


# =======================================================



# test spectrogram
# with real input matlab outputs a 1-sided PSD
x0 = readdlm(joinpath("../test/data", "spectrogram_x.txt"),'\t')
f0 = readdlm(joinpath("../test/data", "spectrogram_f.txt"),'\t')
t0 = readdlm(joinpath("../test/data", "spectrogram_t.txt"),'\t')
p0 = readdlm(joinpath("../test/data", "spectrogram_p.txt"),'\t')
# p, t, f = spectrogram(x0, n=256, r=10, m=128)
pm = mean(p0,2)
P = Periodogramt(x0,n=256,noverlap=128,fs=10,twosided=false,window=false)
Pp = power(P)
Ps, Pt, Pf = stf(P)

@test_approx_eq Pp pm
@test_approx_eq f0 Pf
@test_approx_eq t0 Pt
@test_approx_eq p0 Ps

# # ~~~~ TESTS FOR DSP.Periodogram.welch_pgram ~~~~
data = Float64[0:7]

# ~~~~~~~~~~~ This one tests periodogram ~~~~~~~~~~
#Matlab: p = pwelch(0:7, [1, 1, 1, 1, 1, 1, 1, 1], 0, 8, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, length(data), 0)) Float64[ 98.0,
                                                             13.656854249492380,
                                                              4.0,
                                                              2.343145750507620,
                                                              2.0,
                                                              2.343145750507620,
                                                              4.0,
                                                             13.656854249492380]

# ~~~~~~~~ Tests with no window ~~~~~~~~~~~~~~~~~~~
# Matlab: p = pwelch(0:7, [1, 1], 0, 2, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 2, 0)) Float64[34.5, 0.5]

# Matlab: p = pwelch(0:7, [1, 1, 1], 0, 3, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 3, 0)) Float64[25.5, 1.0, 1.0]

# Matlab: p = pwelch(0:7, [1, 1, 1], 1, 3, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 3, 1)) Float64[35.0, 1.0, 1.0]

# Matlab: p = pwelch(0:7, [1, 1, 1, 1], 1, 4, 1, 'twosided')
@test_approx_eq power(welch_pgram(data, 4, 1)) Float64[45, 2, 1, 2]

# ~~~~~~~~~ Tests with window ~~~~~~~~~~~~~~~
data = Float64[0:7]

# ~~~~~~~~~~~ This one tests periodogram ~~~~~~~~~~~~
# ~ If functionality of the other arguments has been
# ~ tested above, we only test here that the correct
# ~ value of the spectral density is obtained when
# ~ using a window. More tests to be added if needed
#Matlab: p = pwelch(0:7, window_func(8), 0, 8, 1, 'twosided')
cases = {hamming => Float64[65.461623986801527,
20.556791795515764,
0.369313143650544,
0.022167446610882,
0.025502985564107,
0.022167446610882,
0.369313143650544,
20.556791795515764],
bartlett => Float64[62.999999999999993,
                            21.981076052592442,
                             0.285714285714286,
                             0.161781090264695,
                             0.142857142857143,
                             0.161781090264695,
                             0.285714285714286,
                            21.981076052592442]	}

for (window, expected) in cases
    @test_approx_eq power(welch_pgram(data, length(data), 0, window)) expected
end

end


