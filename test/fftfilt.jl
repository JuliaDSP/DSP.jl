using DSP, Base.Test

for xlen in 2.^(7:18).-1, blen in 2.^(1:6).-1
	b = randn(blen)
	x = rand(xlen)
	filtres = filt(b, [1.0], x)
	fftres = fftfilt(b, x)
	firres  = firfilt(b, x)
	@test_approx_eq filtres fftres
	@test_approx_eq filtres firres
end
