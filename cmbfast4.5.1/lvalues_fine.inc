c have l's in l(j); from 2-10,12,15,20,30,40,60,80,100, every 50 beyond
c This are l's which will be calculated, the rest will be interpolated


c WARNING: If you want to use an even finer sampling, be careful that
c the arrays are correctly dimensioned. You must make sure that the
c number of ls in the array is lower than, lmax=20+l0max/10 . If this
c is not the case you must change the formula for lmax so that it is
c sufficiently large. lmax appears in several places.

        lind=1
        do 22 lvar=2,10
           l(lind)=lvar
           lind=lind+1
 22     continue
        l(lind)=12
        lind=lind+1
        l(lind)=15

        do 24 lvar=20,lmo,10
           lind=lind+1
           l(lind)=lvar
 24     continue
        l0=lind
