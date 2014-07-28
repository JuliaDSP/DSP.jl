:mod:`Windows` - window functions
=================================

.. function:: rect(n)

    Rectangular window function of length ``n``.

.. function:: hanning(n)

    Hanning window of length ``n``.

.. function:: hamming(n) 

    Hamming window of length ``n``.

.. function:: tukey(n, alpha)

    Tukey window of length ``n``, parameterized by ``alpha``. For
    ``alpha`` = 0, the window is equivalent to a rectangular window.
    For ``alpha`` = 1, the window is a Hann window.

.. function:: cosine(n)

    Cosine window of length ``n``. Also called the sine window for
    obvious reasons.

.. function:: lanczos(n)

    Lanczos window of length ``n``.

.. function:: triang(n)

    Triangular window of length ``n``.

.. function:: bartlett(n)

    Bartlett window of length ``n``.

.. function:: gaussian(n, sigma)

    Gaussian window of length ``n`` parameterized by the standard
    deviation ``sigma``.

.. function:: bartlett_hann(n)

    Bartlett-Hann window of length ``n``.

.. function:: blackman(n)

    "Exact" Blackman window, alpha = 0.16.

.. function:: kaiser(n, alpha)

    Kaiser window of length ``n`` parameterized by ``alpha``.

.. function:: dpss(n, nw, ntapers=iceil(2*nw)-1)

    The first ``ntapers`` discrete prolate spheroid sequences (Slepian
    tapers) as an ``n`` × ``ntapers`` matrix. The signs of the tapers
    follow the convention that the first element of the skew-symmetric
    (odd) tapers is positive. The time-bandwidth product is given by
    ``nw``.

.. function:: dpsseig(A, nw)

    Eigenvalues of the DPSS matrix, representing the ratios of the
    power within the main lobe to power in the sidelobe (i.e. leakage).
    ``A`` is the output of :func:`dpss`, and ``nw`` is the
    time-bandwidth product provided to :func:`dpss` as input.
