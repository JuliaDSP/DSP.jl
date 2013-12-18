module Util

export unwrap!, unwrap

function unwrap!(p)
    if length(p) < 2
        return p
    end
    for i = 2:length(p)
        d = p[i] - p[i-1]
        if abs(d) > pi
            p[i] -= floor((d+pi) / (2pi)) * 2pi
        end
    end
    return p
end

unwrap(p) = unwrap!(copy(p))

end # end module definition
