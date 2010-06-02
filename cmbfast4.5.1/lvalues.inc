c have l's in l(j); from 2-10,12,15,20,30,40,60,80,100, every 50 beyond
c This are l's which will be calculated, the rest will be interpolated
        lind=1
        do 22 lvar=2,10
           l(lind)=lvar
           lind=lind+1
 22     continue
        l(lind)=12
        lind=lind+1
        l(lind)=15
        lind=lind+1
        l(lind)=20
        lind=lind+1
        l(lind)=30
        lind=lind+1
        l(lind)=40
        lind=lind+1
        l(lind)=50
        lind=lind+1
        l(lind)=60
        lind=lind+1
        l(lind)=70
        lind=lind+1
        l(lind)=80
        lind=lind+1
        l(lind)=90
        lind=lind+1
        l(lind)=110
        lind=lind+1
        l(lind)=130

        do 24 lvar=150,lmo,50
           lind=lind+1
           l(lind)=lvar
 24     continue
        l0=lind
