c
c       MPI COMMUNICATION SUBROUTINES
c
c       These subroutines were written by David Spergel
c
c
c
c

#ifdef MPIBUILD

      SUBROUTINE TRANSFER_CL(cl,cpl,ccl,ckkl,ctkl,ctl,ctel,ctbl,ctcl,
     >                  itflag,nn,l0)
c
c       this routine accumulates the values of the cl's from all of the slave
c       processors
c
      include 'cmbfast.inc'
      INCLUDE 'mpif.h'
      integer myid,sender,numprocs,ierr,status(MPI_STATUS_SIZE)

      INTEGER itflag,nn,l0
      DOUBLE PRECISION cl(l0,nn),cpl(l0,nn),ccl(l0,nn),ckkl(l0,nn),
     &                 ctkl(l0,nn)
      DOUBLE PRECISION ctl(l0,nn),ctel(l0,nn),ctbl(l0,nn),ctcl(l0,nn)
c
      INTEGER in,j,icount,iproc
      DOUBLE PRECISION datat(l0max*nnmax*9)
c

      CALL MPI_COMM_RANK_( MPI_COMM_WORLD, myid, ierr )          !MPI
      CALL MPI_COMM_SIZE_( MPI_COMM_WORLD, numprocs, ierr )
c
c

      IF(MYID.ne.0) THEN
        icount = 0
        do in=1,nn
        do j = 1,l0
                if(itflag.ne.2) then
                        datat(icount+j) = cl(j,in)

                        datat(icount+l0+j)= cpl(j,in)
                        datat(icount+2*l0+j) = ccl(j,in)
                        datat(icount+3*l0+j) = ckkl(j,in)
                        datat(icount+4*l0+j) = ctkl(j,in)
                endif
        end do
        if(itflag.ne.2) icount = icount+5*l0
        do j = 1,l0
                if(itflag.ne.0) then
                        datat(icount+j) = ctl(j,in)
                        datat(icount+l0+j)= ctel(j,in)
                        datat(icount+2*l0+j) = ctbl(j,in)
                        datat(icount+3*l0+j) = ctcl(j,in)
                endif
        end do
        if(itflag.ne.0) icount = icount+4*l0
        end do
        CALL MPI_SEND_(datat,icount,MPI_DOUBLE_PRECISION
     >                  ,0,0,MPI_COMM_WORLD,IERR)
      ELSE
        DO iproc = 1,numprocs-1
                icount = 0
                if(itflag.ne.2)  icount = icount+5*l0*nn
                if(itflag.ne.0)  icount = icount+4*l0*nn
                CALL MPI_RECV_(datat,icount,MPI_DOUBLE_PRECISION,
     >                                  MPI_ANY_SOURCE,MPI_ANY_TAG,
     >                                  MPI_COMM_WORLD,STATUS,IERR)
                icount = 0
                do in = 1,nn
                do j = 1,l0
                if(itflag.ne.2) then
                        cl(j,in) = datat(icount+j) + cl(j,in)
                        cpl(j,in) = datat(icount+l0+j)+ cpl(j,in)
                        ccl(j,in) = datat(icount+2*l0+j)+ ccl(j,in)
                        ckkl(j,in) = datat(icount+3*l0+j) + ckkl(j,in)
                        ctkl(j,in) = datat(icount+4*l0+j) + ctkl(j,in)
                endif
                end do
                if(itflag.ne.2) icount = icount+5*l0
                do j = 1,l0
                if(itflag.ne.0) then
                        ctl(j,in) = datat(icount+j) +ctl(j,in)
                        ctel(j,in)=datat(icount+l0+j)+ ctel(j,in)
                        ctbl(j,in)=datat(icount+2*l0+j) + ctbl(j,in)
                        ctcl(j,in)=datat(icount+3*l0+j) + ctcl(j,in)
                endif
                end do
                if(itflag.ne.0) icount = icount+4*l0
                end do
         END DO
        ENDIF
        RETURN
        END
c
c
c
c
      SUBROUTINE COMMUNICATE(d,dp,dkk,dt,dte,dtb,nstep,nk,itflag)
c
c       this shares the values of d,dp,dkk,dt,dte, and dtb
c       with the other processorss
c
      include 'cmbfast.inc'
      INCLUDE 'mpif.h'
      integer myid,sender,numprocs,ierr,status(MPI_STATUS_SIZE)

      integer partner
      external partner
      integer nk,nstep,itflag
      integer istep,iproc
      double precision d(nk0,nstep0),dp(nk0,nstep0),dkk(nk0,nstep0)
      double precision dt(nk0,nstep0),dte(nk0,nstep0),dtb(nk0,nstep0)
      double precision datat(6*nstep0*nk0),datar(6*nstep0*nk0)
c
      integer i,ik,icount,icount0,isender
      integer ncommstep,nbit

      CALL MPI_COMM_RANK_( MPI_COMM_WORLD, myid, ierr )          !MPI
      CALL MPI_COMM_SIZE_( MPI_COMM_WORLD, numprocs, ierr )
c
c       find maximum value of nstep
c
        IF(MYID.eq.0) THEN
                nstepm = nstep
                do iproc = 1,numprocs-1
                        CALL MPI_RECV_(nsteps,1,
     >                                  MPI_INTEGER,
     >                                  MPI_ANY_SOURCE,MPI_ANY_TAG,
     >                                  MPI_COMM_WORLD,STATUS,IERR)
                        nstepm = max(nstepm,nsteps)
                end do
        ELSE
                CALL MPI_SEND_(nstep,1,MPI_INTEGER,0,
     >                                  0,MPI_COMM_WORLD,IERR)
c                write(*,*)'nstep=',nstep
        ENDIF
        CALL MPI_BCAST_(nstepm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

c
c
c       Store data for transfer in 1d file
c
        icount = 0
        do ik=myid+1,nk,numprocs
                if(itflag.ne.2) then
                                do i = 1,nstepm
                                        datat(i+icount) = d(ik,i)
                                        datat(i+icount+nstepm) = dp(ik,i)
                                        datat(i+icount+2*nstepm) = dkk(ik,i)
                                end do
                                icount = icount+3*nstepm
                endif
                if(itflag.ne.0) then
                                do i = 1,nstepm
                                        datat(i+icount) = dt(ik,i)
                                        datat(i+icount+nstepm) = dte(ik,i)
                                        datat(i+icount+2*nstepm) = dtb(ik,i)
                                end do
                                icount = icount+3*nstepm
                endif
        end do
        icount0 = icount

c
c       compute # of bits in binary representation of numprocs
c
        DO nbit = 1,10
                if(2**nbit.ge.numprocs) go to 1990
        END DO
1990    ncommstep = 2**nbit-1

c
c
c       send data to other processors
c
        DO istep = 1,ncommstep
                CALL MPI_BARRIER_(MPI_COMM_WORLD,ierr)
                iproc = partner(MYID,istep,numprocs)
                IF(iproc.le.numprocs-1) THEN

                        IF(iproc.lt.MYID) then

c       write(*,*) 'sending data from',myid,'to',iproc
                                icount =0
                                if(itflag.ne.2) icount = icount+3*nstepm*
     >                                  (nk/numprocs+1)
                                if(itflag.ne.0) icount = icount+3*nstepm*
     >                                  (nk/numprocs+1)
                                CALL MPI_SEND_(datat,icount,
     >                                  MPI_DOUBLE_PRECISION,iproc,
     >                                  0,MPI_COMM_WORLD,IERR)
                                CALL MPI_RECV_(datar,icount,
     >                                  MPI_DOUBLE_PRECISION,
     >                                  MPI_ANY_SOURCE,MPI_ANY_TAG,
     >                                  MPI_COMM_WORLD,STATUS,IERR)
                                isender = STATUS(MPI_SOURCE)
                                icount = 0
                                do ik = isender+1,nk,numprocs
                                if(itflag.ne.2) then
                                do i = 1,nstepm
                                        d(ik,i) = datar(i+icount)
                                        dp(ik,i) = datar(i+icount+nstepm)
                                        dkk(ik,i) = datar(i+icount+2*nstepm)
                                end do
                                icount = icount+3*nstepm
                                endif
                                if(itflag.ne.0) then
                                do i = 1,nstepm
                                        dt(ik,i) = datar(i+icount)
                                        dte(ik,i) = datar(i+icount+nstepm)
                                        dtb(ik,i) = datar(i+icount+2*nstepm)
                                end do
                                icount = icount+3*nstepm
                                endif
                                end do
                        ELSE
                                icount =0
                                if(itflag.ne.2) icount = icount+3*nstepm*
     >                                  (nk/numprocs+1)
                                if(itflag.ne.0) icount = icount+3*nstepm*
     >                                  (nk/numprocs+1)
                                CALL MPI_RECV_(datar,icount,
     >                                  MPI_DOUBLE_PRECISION,
     >                                  MPI_ANY_SOURCE,MPI_ANY_TAG,
     >                                  MPI_COMM_WORLD,STATUS,IERR)
                                isender = STATUS(MPI_SOURCE)

                                icount = 0
                                do ik = isender+1,nk,numprocs
                                if(itflag.ne.2) then
                                do i = 1,nstepm
                                        d(ik,i) = datar(i+icount)
                                        dp(ik,i) = datar(i+icount+nstepm)
                                        dkk(ik,i) = datar(i+icount+2*nstepm)
                                end do
                                icount = icount+3*nstepm
                                endif
                                if(itflag.ne.0) then
                                do i = 1,nstepm
                                        dt(ik,i) = datar(i+icount)
                                        dte(ik,i) = datar(i+icount+nstepm)
                                        dtb(ik,i) = datar(i+icount+2*nstepm)
                                end do
                                icount = icount+3*nstepm
                                endif
                                end do
c                       write(*,*) 'sending data from',myid,'to',iproc
                                icount =0
                                if(itflag.ne.2) icount = icount+3*nstepm*
     >                                  (nk/numprocs+1)
                                if(itflag.ne.0) icount = icount+3*nstepm*
     >                                  (nk/numprocs+1)
                                CALL MPI_SEND_(datat,icount,
     >                                  MPI_DOUBLE_PRECISION,iproc,
     >                                  0,MPI_COMM_WORLD,IERR)
                        ENDIF
                ENDIF
        END DO
        RETURN
        END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER FUNCTION partner(myid,istep,numprocs)
      INTEGER, INTENT(IN) :: myid,istep,numprocs
      INTEGER ibit
      partner =0
      do ibit = 1,10
        if(mod((myid)/2**(ibit-1),2).ne.mod(istep/2**(ibit-1),2))
     >          partner=partner+2**(ibit-1)
        if(2**ibit.ge.numprocs) return
      end do
      RETURN
      END

#endif
