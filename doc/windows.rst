:mod:`Windows` - window functions
=================================

.. function:: rect(n::Integer)

    Rectangular window function of length n.

.. function:: hanning(n::Integer)

    Hanning window of length n.

.. function:: hamming(n::Integer) 

    Hamming window of length n.

.. function:: tukey(n::Integer, alpha::Real)

    Tukey window of length n, parameterized by alpha. For alpha = 0, the window is equivalent to a rectangular window. For alpha = 1, the window is a Hann window.

.. function:: cosine(n::Integer)

    Cosine window of length N. Also called the sine window for obvious reasons.

.. function:: lanczos(n::Integer)

    Lanczos window of length n.

.. function:: triang(n::Integer)

    Triangular window of length n.

.. function:: bartlett(n::Integer)

    Bartlett window of length n.

.. function:: gaussian(n::Integer, sigma::Real)

    Gaussian window of length N parameterized by the standard deviation sigma

.. function:: bartlett_hann(n::Integer)

    Bartlett-Hann window of length n.

.. function:: blackman(n::Integer)

    "Exact" Blackman window, alpha=0.16.

.. function:: kaiser(n::Integer, alpha::Real)

    Kaiser window parameterized by alpha.

.. function:: dpss(n::Int, nw::Real, ntapers::Int=iceil(2*nw)-1)

    The first ``ntapers`` discrete prolate spheroid sequences (Slepian tapers) as an ``n`` x ``ntapers`` matrix. The signs of the tapers follow the convention that the first element of the skew-symmetric (odd) tapers is positive. The time-bandwidth product is given by ``nw``.
