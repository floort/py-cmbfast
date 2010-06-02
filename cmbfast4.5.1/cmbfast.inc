c     l0max is used to dimension arrays. It must be bigger than the
c     maximum l to be calculated by at least 300.
c
        parameter (l0max=5300)

c
c       Dimension of the l arrays.
c WARNING: If you want to use a finer sampling in l than in
c lvalues_fine.inc, you must make sure that the
c number of ls in the array is lower than lmax=20+l0max/10 . If this
c is not the case you must change the formula for lmax so that it
c becomes large enough.

        parameter (lmax=20+l0max/10)


c       nnmax sets the maximum number of spectral indeces that
c       can be calculated.

        parameter(nnmax=11)

c       maximum number of elements in (a,w) table for dark energy models
        parameter(nq0=2000)
c       ntfmax sets the maximum number of transfer functions that
c       can be calculated.

        parameter(ntfmax=211)
        integer nthermo
        parameter(nthermo=10000)

c ntau sets the time sampling for lensing
        parameter (ntau=201)

c
c       nk0 and nstep0 set the dimension of the arrays
c       of k and time values were the sources are calculated

        parameter (nk0=300)
c       parameter (nstep0=2400)
        parameter (nstep0=7000)

