module Extensions
export extend, extend!,
    ZeroExt, ConstantExt, LinearExt, PeriodicExt

abstract AbstractExtend

# ( ... x[n] | 0 0 ... )
immutable ZeroExtend        <: AbstractExtend end # fill with 0
const ZeroExt = ZeroExtend()

# ( ... x[n] | x[n] x[n] ... )
immutable ConstantExtend    <: AbstractExtend end # 0th-order extrapolation
const ConstantExt = ConstantExtend()

# ( ... x[n-1] x[n] | 2*x[n]-x[n-1] ... )
immutable LinearExtend      <: AbstractExtend end # 1st-order extrapolation
const LinearExt = LinearExtend()
 
# ( ... x[n] | x[1] x[2] ... )
immutable PeriodicExtend    <: AbstractExtend end
const PeriodicExt = PeriodicExtend()

# ( ... x[n] | x[n] x[n-1] ... )
immutable SymmetricExtend   <: AbstractExtend end # mirror extension
const SymmetricExt = SymmetricExtend()

# ( ... x[n] | x[n-1] x[n-2] ... )
immutable Symmetric0Extend  <: AbstractExtend end # without duplicates at boundary
const Symmetric0Ext = SymmetricExtend()


# copies s to out[offset+1:offset+length(s)]
function extend{T}(n::Int, s::AbstractVector{T}, args...)
    out = Array(T, n)
    extend!(out, s, args...)
end
function extend!{T}(out::AbstractVector{T}, s::AbstractVector{T}, ext::AbstractExtend, offset::Integer=0)
    ns = length(s)
    nout = length(out)
    0 <= offset || throw(ArgumentError("negative offset"))
    ns+offset <= nout || throw(BoundsError())
    
    copy!(out, offset+1, s, 1, ns)
    addextensions!(out, ns, ext, offset)
    return out
end
function extend!{T}(out::AbstractVector{T}, s::AbstractVector{T}, ext::AbstractExtend, se::Function, offset::Integer=0)
    extend!(out, s, ext, offset)
    smoothedges!(out, length(s), se, offset)
    return out
end
# apply window function to extended ranges
function smoothedges!{T}(out::AbstractVector{T}, ns::Integer, se::Function, offset::Integer)
    nout = length(out)
    
    win = se(offset<<1)
    for i = 1:offset
        @inbounds out[i] *= win[i]
    end
    nsoff = -ns-offset
    win = se((nout+nsoff)<<1)
    nsoff = (nsoff<<1)-1+nout
    for i = ns+offset+1:nout
        @inbounds out[i] *= win[i+nsoff]
    end
    return out
end

## EXTENSION IMPLEMENTATIONS
# unsafe, no checks for ns or offset

function addextensions!{T}(out::AbstractVector{T}, ns::Integer, ::ZeroExtend, offset::Integer)
    nout = length(out)
    
    z = zero(T)
    for i = 1:offset
        @inbounds out[i] = z
    end
    for i = ns+offset+1:nout
        @inbounds out[i] = z
    end
    return out
end
function addextensions!{T}(out::AbstractVector{T}, ns::Integer, ::ConstantExtend, offset::Integer)
    nout = length(out)
    ns > 0 || return addextensions!(out, ns, ZeroExt, offset)
    
    a = out[offset+1]
    for i = 1:offset
        @inbounds out[i] = a
    end
    a = out[ns+offset]
    for i = ns+offset+1:nout
        @inbounds out[i] = a
    end
    return out
end
function addextensions!{T}(out::AbstractVector{T}, ns::Integer, ::LinearExtend, offset::Integer)
    nout = length(out)
    ns > 1 || return addextensions!(out, ns, ConstantExt, offset)
    
    k = 1+offset
    val = out[k]
    d = val - out[k+1]
    for i = offset:-1:1
        @inbounds out[i] = val + d*(k-i)
    end
    
    k = ns+offset
    val = out[k]
    d = val - out[k-1]
    for i = ns+offset+1:nout
        @inbounds out[i] = val + d*(i-k)
    end
    return out
end
function addextensions!{T}(out::AbstractVector{T}, ns::Integer, ::PeriodicExtend, offset::Integer)
    nout = length(out)
    ns > 0 || return addextensions!(out, ns, ZeroExt, offset)
    
    offsetp1 = offset + 1   
    if offset <= ns  # more likely -> faster loop
        for i = 1:offset
            @inbounds out[i] = out[i+ns]
        end
    else
        for i = 1:offset
            @inbounds out[i] = out[offsetp1 + mod(i-offsetp1, ns)]
        end
    end
    if nout-(ns+offset) <= ns  # more likely -> faster loop
        for i = ns+offset+1:nout
            @inbounds out[i] = out[i-ns]
        end
    else
        for i = ns+offset+1:nout
            @inbounds out[i] = out[offsetp1 + mod(i-offsetp1, ns)]
        end
    end
    return out
end

end # module

