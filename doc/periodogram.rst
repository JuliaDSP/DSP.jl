:mod:`Periodogram` - periodogram estimation
===========================================

.. function:: arraysplit(s, n::Integer, m::Integer)

    Split an array into arrays of length n with overlapping regions of length m.

.. function:: periodogram(s)

    Compute periodogram of a signal by FFT.

.. function:: welch_pgram(s, n, m)

    Compute Welch periodogram of a signal based on n segments with overlap m.

.. function:: bartlett_pgram(s, n)

    Compute Bartlett periodogram. This is equivalent to welch_pgram(s, n, 0).
