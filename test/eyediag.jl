using Test, DSP

@testset "Eye Diagram Tests" begin #Scope for test data

testbits = [0,1,0,1,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0]
#@show length(testbits)

#Simple function to build a triangular bit pattern from a bit sequence
function BuildTriangPat(bitseq; tbit::Float64=1.0)
	y = 1.0 .* bitseq #Get floating point values
	x = collect((0:1:(length(bitseq)-1)) .* tbit)
	return (x, y)
end

#For debug purposes:
function dbg_showeyedata(eyedata)
	for i in 1:length(eyedata.vx)
		println("$i:")
		@show eyedata.vx[i]
		@show eyedata.vy[i]
	end
end

tbit = 1.0
(x,y) = BuildTriangPat(testbits, tbit = tbit)

@testset "Eye Diagram: Centered" begin
	eyedata = buildeye(x, y, tbit, 2.0*tbit, tstart=0.0*tbit)

	@test length(eyedata.vx) == length(eyedata.vy)
	@test length(eyedata.vx) == length(testbits) - 1

#	dbg_showeyedata(eyedata)

	for i in 1:(length(eyedata.vx)-1)
		@test eyedata.vx[i] == [0.0, 1.0, 2.0]
		@test eyedata.vy[i] == Float64[testbits[i], testbits[i+1], testbits[i+2]]
	end

	i = length(eyedata.vx)
	@test eyedata.vx[end] == [0.0, 1.0] #Cannot extrapolate past last data point
	@test eyedata.vy[i] == Float64[testbits[i], testbits[i+1]]
end

@testset "Eye Diagram: Early" begin
	eyedata = buildeye(x, y, tbit, 2.0*tbit, tstart=0.8*tbit)

	@test length(eyedata.vx) == length(eyedata.vy)
	@test length(eyedata.vx) == length(testbits) - 2 #Skipped 1st data point

#	dbg_showeyedata(eyedata)

	for i in 1:length(eyedata.vx)
		@test eyedata.vx[i] ≈ [0.2, 1.2]
		@test eyedata.vy[i] == Float64[testbits[i+1], testbits[i+2]]
	end
end

@testset "Eye Diagram: Late" begin
	eyedata = buildeye(x, y, tbit, 2.0*tbit, tstart=0.2*tbit)

	@test length(eyedata.vx) == length(eyedata.vy)
	@test length(eyedata.vx) == length(testbits) - 2 #Skipped 1st data point

#	dbg_showeyedata(eyedata)

	for i in 1:length(eyedata.vx)
		#Last data point in sweep should be >= 2.0*tbit:
		@test eyedata.vx[i] ≈ [0.8, 1.8]
		@test eyedata.vy[i] == Float64[testbits[i+1], testbits[i+2]]
	end
end
end
