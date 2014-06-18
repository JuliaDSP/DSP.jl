
module Periodograms

export Periodogramt, power, freq, spectrogram, time, periodogram, welch_pgram, bartlett_pgram

typealias DSPNumber FFTW.fftwNumber
typealias DSPReal FFTW.fftwReal
typealias DSPComplex FFTW.fftwComplex

type Periodogramt{T<:DSPNumber}   # Periodogram is the module name, so there is a conflict
    s::Array{T,1}
    n::Integer                      # 0 < n <= length(s)
    m::Integer                      # 0 <= m < n
    r::Real                         # positive
    window::Union(Function,Bool)    # false || Function
    sided::Integer                  # 1 || 2
end
function Periodogramt{T<:DSPNumber}(s::Array{T,1}; n = length(s), m = 0, r = 1, window = false, sided = 1)
    po = Periodogramt(s,n,m,r,window,sided)
    errorcheck(po)
    return po
end
function Periodogramt{T<:DSPNumber}(s::Array{T,2}; args...)
    if size(s,1)==1 || size(s,2)==1
        return Periodogramt(vec(s);args...)  # convert to a vector
    else
        error("not supported for 2D arrays")
    end
end
function errorcheck{T<:DSPNumber}(P::Periodogramt{T})
    !(0 < P.n <= length(P.s)) && error("use 0 < n <= length(s)")
    !(0 <= P.m < P.n) && error("use 0 <= m < n")
    P.r <= 0 && error("r is not positive")
    P.window == true && error("window = true not supported")
    P.sided != 1 && P.sided != 2 && error("sided = ",P.sided," not supported")
end

# ======= the old DSP functions =======

# sided=2 default for now since that is the current implementation in DSP
# could be changed to call directly power(Periodogramt...)
function periodogram{T<:DSPNumber}(s::Array{T,1})
    Periodogramt(s,sided=2)
end
function periodogram{T<:DSPNumber}(s::Array{T,1}, window)
    Periodogramt(s,sided=2,window=window)
end
function welch_pgram{T<:DSPNumber}(s::Array{T,1}, n, m )
    Periodogramt(s,n=n,m=m,sided=2)
end
function welch_pgram{T<:DSPNumber}(s::Array{T,1}, n, m, window )
    Periodogramt(s,n=n,m=m,window=window,sided=2)
end
function bartlett_pgram{T<:DSPNumber}(s::Array{T,1}, n)
    welch_pgram(s, n, 0)
end
function bartlett_pgram{T<:DSPNumber}(s::Array{T,1}, n, window)
    welch_pgram(s, n, 0, window)
end
function spectrogram{T<:DSPNumber}(s::Array{T,1}; n=int(length(s)/8), m=int(n/2), r=1, w=(n)->ones(n))
    Periodogramt(s,n=n,m=m,r=r,window=w,sided=2)
end

# ======= methods to apply to type object =======

# spectrum vector
function power{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    sided = P.sided
    eltype(P.s)<:DSPComplex && (sided = 2)
    if length(P.s) == P.n
        return power_p(P.s,P.r,P.window,sided)
    else
        return power_w(P.s,P.n,P.m,P.r,P.window,sided)
    end
end
# frequency vector
function freq{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    sided = P.sided
    eltype(P.s)<:DSPComplex && (sided = 2)
    return periodogramfreq(P.n, r=P.r, sided = sided)
end
# spectrogram matrix
function spectrogram{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    sided = P.sided
    eltype(P.s)<:DSPComplex && (sided = 2)
    return sp_gram(P.s,P.n,P.m,P.r,P.window,sided)
end
# time vector (for spectrogram)
function time{T<:DSPNumber}(P::Periodogramt{T})
    errorcheck(P)
    index = arrayspliti(P.s, P.n, P.m)
    return ( (0:length(index)-1)*(P.n-P.m) + P.n/2) / P.r
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

# does not call power_p in the loop to avoid tmp array allocation
# welch type spectrum
function power_w{T<:DSPNumber}(s::Array{T,1}, n::Integer, m::Integer, r::Real, window::Union(Function,Bool), sided::Integer)
    tsReal = typeof(real(s[1]))
    tsComplex = typeof(complex(s[1]))
    if isa(window,Function)
        w = window(n)
    end
    index = arrayspliti(s, n, m)
    if sided == 2
        p = zeros(tsReal,n)
        split = Array(tsComplex,n)
        for i in index
            if isa(window,Function)
                split[:] = s[i:i+n-1].*w
            else
                split[:] = s[i:i+n-1]
            end
            periodogramtl!(split)
            p += real(split)
        end
    else #sided == 1
        np = div(n,2)+1
        p = zeros(tsReal,np)
        split = Array(tsReal,n)
        splitrf = Array(tsComplex, np)
        for i in index
            if isa(window,Function)
                split[:] = s[i:i+n-1].*w
            else
                split[:] = s[i:i+n-1]
            end
            periodogramtl1s!(splitrf, split)
            p += real(splitrf)
        end
    end
    if isa(window,Function)
        p /= (r*length(index)*sum(w.^2))
    else
        p /= (r*length(index)*n)
    end
    return p
end

# periodogram spectrum
function power_p{T<:DSPNumber}(s::Array{T,1}, r::Real, window::Union(Function,Bool), sided::Integer)
    tsReal = typeof(real(s[1]))
    tsComplex = typeof(complex(s[1]))
    n=length(s)
    if isa(window,Function)
        w = window(n)
    end
    if sided == 2
        pc = Array(tsComplex, n)
        p = Array(tsReal, n)
        if isa(window,Function)
            pc[:] = s.*w
        else
            pc[:] = s
        end
        periodogramtl!(pc)
    else #sided == 1
        np = div(n,2)+1
        #p = Array(typeof(s[1]*im), np)
        pc = Array(tsComplex, np)
        p = Array(tsReal, np)
        if isa(window,Function)
            periodogramtl1s!(pc, s.*w)
        else
            periodogramtl1s!(pc, s)
        end
    end
    if isa(window,Function)
        p[:] = pc/(r*sum(w.^2))
    else
        p[:] = pc/(r*n)
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
function sp_gram{T<:DSPNumber}(s::Array{T,1}, n::Integer, m::Integer, r::Real, window::Union(Function,Bool), sided::Integer)
    tsReal = typeof(real(s[1]))
    tsComplex = typeof(complex(s[1]))
    if isa(window,Function)
        w = window(n)
    end
    index = arrayspliti(s, n, m)
    if sided == 2
        p = Array(tsReal,n,length(index))
        split = Array(tsComplex,n)
        j = 1
        for i in index
            if isa(window,Function)
                split[:] = s[i:i+n-1].*w
            else
                split[:] = s[i:i+n-1]
            end
            periodogramtl!(split)
            p[:,j] = real(split)
            j += 1
        end
    else #sided == 1
        np = div(n,2)+1
        p = zeros(tsReal,np,length(index))
        split = Array(tsReal,n)
        splitrf = Array(tsComplex, np)
        j = 1
        for i in index
            if isa(window,Function)
                split[:] = s[i:i+n-1].*w
            else
                split[:] = s[i:i+n-1]
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

# frequency vector from window length n and sampling rate r, sided is either 1 (default) or 2
# r = n gives a vector of wavenumbers
function periodogramfreq{T<:Integer}(n::T; r::Real = 1, sided::T = 1)
    if sided == 1
        nf = div(n,2)+1  # the highest frequency + 1
    elseif sided == 2
        nf = n
    else
        error("sided = ", sided, " not supported")
    end
    f = Array(typeof(1.0),nf)
    rn=r/n
    for i = 1:length(f)
        if i <= div(n,2)+1
            @inbounds f[i] = (i-1)*rn
        else
            @inbounds f[i] = (-N+i-1)*rn
        end
    end
    return f
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
# initialize variables before timing
P=Periodogramt(x0,n=nn,m=mm,sided=2,window=n->ones(n))
Pp = power(P)
P.sided=1
Pp1 = power(P)
P.sided=2
wp = Periodogram.welch_pgram(x0,nn,mm)

tn=10

println("  P 2 sided")
@time begin for i=1:tn
Pp = power(P)
end
end
println(typeof(Pp)," ",length(Pp))

println("  old welch")
@time begin for i=1:tn
wp = Periodogram.welch_pgram(x0,nn,mm)
end
end
println(typeof(wp)," ",length(wp))

println("  P 1 sided")
P.sided=1
@time begin for i=1:tn
Pp1 = power(P)
end
end
println(typeof(Pp1)," ",length(Pp1))

@test_approx_eq Pp wp

println("\n periodograms : \n")
nn=N
mm=0
P=Periodogramt(x0,n=nn,m=mm,sided=2)
Pp = power(P)
P.sided=1
Pp1 = power(P)
P.sided=2
wp = Periodogram.periodogram(x0)

println("  P 2 sided")
@time begin for i=1:tn
Pp = power(P)
end
end
println(typeof(Pp)," ",length(Pp))

println("  old periodogram")
@time begin for i=1:tn
wp = Periodogram.periodogram(x0)
end
end
println(typeof(wp)," ",length(wp))

println("  P 1 sided")
P.sided=1
@time begin for i=1:tn
Pp1 = power(P)
end
end
println(typeof(Pp1)," ",length(Pp1))

@test_approx_eq Pp wp


println("\n spectrograms : \n")

Ps = spectrogram(P)
sp,t,f = Periodogram.spectrogram(x0,n=nn,m=mm)

println("  P spectrogram")
P.sided=2
@time begin for i=1:tn
Ps = spectrogram(P)
end
end
println(typeof(Ps)," ",length(Ps))

println("  old spectrogram")
@time begin for i=1:tn
sp,t,f = Periodogram.spectrogram(x0,n=nn,m=mm)
end
end
println(typeof(sp)," ",length(sp))

@test_approx_eq Ps sp

println("  P spectrogram 1 sided")
P.sided=1
@time begin for i=1:tn
Ps = spectrogram(P)
end
end
println(typeof(Ps)," ",length(Ps))




# =======================================================



# test spectrogram
# with real input matlab outputs a 1-sided PSD
x0 = readdlm(joinpath("../test/data", "spectrogram_x.txt"),'\t')
f0 = readdlm(joinpath("../test/data", "spectrogram_f.txt"),'\t')
t0 = readdlm(joinpath("../test/data", "spectrogram_t.txt"),'\t')
p0 = readdlm(joinpath("../test/data", "spectrogram_p.txt"),'\t')
# p, t, f = spectrogram(x0, n=256, r=10, m=128)
pm = mean(p0,2)
P = Periodogramt(x0,n=256,m=128,r=10,sided=1)
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


