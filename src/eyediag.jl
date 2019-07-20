#Eye diagram generation
#-------------------------------------------------------------------------------

export buildeye

#==Types
===============================================================================#
#=TODO: Restructure DataEye to be a vector of (x,y) vector pairs.
        - Less likely to be in an invalid state)
        - But less compatible with most plotting tool API (which tend to accept
          vectors of x & y vectors).
=#     
mutable struct DataEye
	vx::Vector{Vector{Float64}}
	vy::Vector{Vector{Float64}}
end
DataEye() = DataEye([], [])


#==
===============================================================================#
#=TODO: - kwargs tbit, teye?:
          Preferable to use named arguments, but can it be done cleanly without
          being confusing to user?
        - Define algorithm to detect first crossing and center eye
          (select `tstart`) accordingly.
=#

"""
    buildeye(x::Vector, y::Vector, tbit::Number, teye::Number; tstart::Number=0)

Segments the `(x,y)` vector into multiple traces of the eye diagram.

Trace data is stored in a `DataEye` object as vectors of `(x,y)` sub-vectors.

Inputs:
  - tbit: Bit period (s).  Defines "trigger point" of each bit start.
  - teye: Eye period (s).  Defines how much of eye data (along x-axis) algorithm is to collect after "trigger point".  Typically, this is set to `2.0*tbit`.
  - tstart: Time of first "trigger point" (s).

Example plotting with Plots.jl:

    #Assumption: (x, y) data generated here.
    tbit = 1e-9 #Assume data bit period is 1ns.

    #Build eye & use tstart to center data.
    eye = buildeye(x, y, tbit, 2.0*tbit, tstart=0.2*tbit)

    plot(eye.vx, eye.vy)
"""
function buildeye(x::Vector, y::Vector, tbit::Number, teye::Number; tstart::Number=0)
	eye = DataEye()

	i = 1
	#skip initial data:
	while i <= length(x) && x[i] < tstart
		i+=1
	end

	wndnum = 0
	inext = i
	while true
		if inext > length(x); break; end #Nothing else to add
		wndstart = tstart+wndnum*tbit
		nexteye = wndstart+tbit
		wndend = wndstart+teye
		istart = inext
		i = istart
		while i <= length(x) && x[i] < nexteye
			i+=1
		end
		inext = i
		while i <= length(x) && x[i] <= wndend
			i+=1
		end
		if i > length(x)
			i = length(x)
		end
		if x[i] > wndend
			i -= 1
		end
		if i > istart
			push!(eye.vx, x[istart:i].-wndstart)
			push!(eye.vy, y[istart:i])
		end
		wndnum += 1
	end
	return eye
end

#Last line
