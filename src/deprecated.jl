# deprecations in 0.8
import .Convolutions.nextfastfft
@deprecate nextfastfft(ns...) nextfastfft.(ns) false

import .Convolutions.conv
@deprecate (conv(u::AbstractVector{T}, v::AbstractVector{T}, A::AbstractMatrix{T}) where T) conv(u, transpose(v), A)

@deprecate(
    filt!(out, f::PolynomialRatio{:z}, x::AbstractVector, si::AbstractVector),
    filt!(out, DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt!(out, f::PolynomialRatio{:z}, x::AbstractArray{<:Any, N}, si::AbstractArray{<:Any, N}) where {N},
    filt!(out, DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt!(out, f::PolynomialRatio{:z}, x::AbstractArray, si::AbstractVector),
    filt!(out, DF2TFilter(f, repeat(si; outer=(1, size(x)[2:end]...))), x)
)
@deprecate(
    filt(f::PolynomialRatio{:z}, x::AbstractVector, si::AbstractVector),
    filt(DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt(f::PolynomialRatio{:z}, x::AbstractArray{<:Any, N}, si::AbstractArray{<:Any, N}) where {N},
    filt(DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt(f::PolynomialRatio{:z}, x::AbstractArray, si::AbstractVector),
    filt(DF2TFilter(f, repeat(si; outer=(1, size(x)[2:end]...))), x)
)
@deprecate(
    filt!(out, f::SecondOrderSections{:z}, x::AbstractVector, si::AbstractVecOrMat),
    filt!(out, DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt!(out, f::SecondOrderSections{:z}, x::AbstractArray, si::AbstractArray),
    filt!(out, DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt!(out, f::SecondOrderSections{:z}, x::AbstractArray, si::AbstractVecOrMat),
    filt!(out, DF2TFilter(f, repeat(si; outer=(1, 1, size(x)[2:end]...))), x)
)
@deprecate(
    filt(f::SecondOrderSections{:z}, x::AbstractVector, si::AbstractVecOrMat),
    filt(DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt(f::SecondOrderSections{:z}, x::AbstractArray, si::AbstractArray),
    filt(DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt(f::SecondOrderSections{:z}, x::AbstractArray, si::AbstractVecOrMat),
    filt(DF2TFilter(f, repeat(si; outer=(1, 1, size(x)[2:end]...))), x)
)
@deprecate(
    filt!(out, f::Biquad{:z}, x::AbstractVector, si::AbstractVector),
    filt!(out, DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt!(out, f::Biquad{:z}, x::AbstractArray{<:Any, N}, si::AbstractArray{<:Any, N}) where {N},
    filt!(out, DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt!(out, f::Biquad{:z}, x::AbstractArray, si::AbstractVector),
    filt!(out, DF2TFilter(f, repeat(si; outer=(1, size(x)[2:end]...))), x)
)
@deprecate(
    filt(f::Biquad{:z}, x::AbstractVector, si::AbstractVector),
    filt(DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt(f::Biquad{:z}, x::AbstractArray{<:Any, N}, si::AbstractArray{<:Any, N}) where {N},
    filt(DF2TFilter(f, copy(si)), x)
)
@deprecate(
    filt(f::Biquad{:z}, x::AbstractArray, si::AbstractVector),
    filt(DF2TFilter(f, repeat(si; outer=(1, size(x)[2:end]...))), x)
)
@deprecate(
    filt!(out, b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number}, x::AbstractVector, si::AbstractVector),
    filt!(out, DF2TFilter(PolynomialRatio(b, a), copy(si)), x)
)
@deprecate(
    filt!(out, b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number}, x::AbstractArray{<:Any, N}, si::AbstractArray{<:Any, N}) where {N},
    filt!(out, DF2TFilter(PolynomialRatio(b, a), copy(si)), x)
)
@deprecate(
    filt!(out, b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number}, x::AbstractArray, si::AbstractVector),
    filt!(out, DF2TFilter(PolynomialRatio(b, a), repeat(si; outer=(1, size(x)[2:end]...))), x)
)
@deprecate(
    filt(b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number}, x::AbstractVector, si::AbstractVector),
    filt(DF2TFilter(PolynomialRatio(b, a), copy(si)), x)
)
@deprecate(
    filt(b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number}, x::AbstractArray{<:Any, N}, si::AbstractArray{<:Any, N}) where {N},
    filt(DF2TFilter(PolynomialRatio(b, a), copy(si)), x)
)
@deprecate(
    filt(b::Union{AbstractVector, Number}, a::Union{AbstractVector, Number}, x::AbstractArray, si::AbstractVector),
    filt(DF2TFilter(PolynomialRatio(b, a), repeat(si; outer=(1, size(x)[2:end]...))), x)
)
