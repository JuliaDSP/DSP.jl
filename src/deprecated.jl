# deprecations in 0.8
import .Util.nextfastfft
@deprecate nextfastfft(ns...) nextfastfft.(ns) false

@deprecate (conv(u::AbstractVector{T}, v::AbstractVector{T}, A::AbstractMatrix{T}) where T) conv(u, transpose(v), A)
