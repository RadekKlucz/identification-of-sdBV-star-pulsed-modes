C Liczy FT od 0 do 50 c/d z rozdzielczoscia 10 pktow na pik, FT zapisze do
c ft50.trf
c
c Uzycie:
c         jkdft inputfile
c
C TAKES POWER SPECTRUM OF REAL TIME SERIES DATA
C IN FORM (TIME,AMPLITUDE)
C
      real, allocatable, dimension(:) :: t, a, freq, pow, pha
      double precision tin,t0,ttemp,fnyqav,pi
      character infile*40, line*256
c

      call getarg(1,infile)
cjk      print *,' input file:   ', infile
      fmax=50.  ! max frequency
      ofac=10.  ! FT resolution 10 points per peak
      write(*,*)
      
      open (1,file=infile,status='old')
      open (2,file='./output/ft50.trf')
c      open (3,file='pspec.peaks')
      pi=3.141592653589793d0
C skip comment lines
      iskip = 0
      do
        read (1, '(A)') line
        if ( line (1:1) /= "#" ) exit
        iskip = iskip + 1
      enddo
      backspace (1)
c set zero time at the first data point; important
c as this is all single precision!
      read (1,*,end=30) t0
C count data lines
      icnt = 0
      do
        read (1,*,end=15) tin
        icnt = icnt + 1
      enddo
15    print *,' Number of entry points: ',icnt
C allocate memory
      allocate (t(icnt), STAT=ierror)
      if ( ierror /= 0 ) then
        write(*,*)"Error trying to allocate t(icnt)"
        stop
      end if
      allocate (a(icnt), STAT=ierror)
      if ( ierror /= 0 ) then
        write(*,*)"Error trying to allocate a(icnt)"
        stop
      end if
      rewind (1)
C skip comment lines
      if ( iskip.gt.0 ) then
        do i=1,iskip
          read (1, '(A)') line
        enddo
      endif
C read in the data
      do i=1,icnt
        read (1,*,end=30) tin,a(i)
        ttemp=tin-t0
        t(i)=ttemp
      enddo
30    np=i-1
c compute useful things
      fnyqav=(np-1)/(2.0*(t(np)-t(1)))
cjk      print '( "Mean Nyquist frequency=",1pe12.3)',fnyqav
cjk      print '( "T(1) =",1pe12.3)',t(1)
cjk      print '( "T(N) =",1pe12.3)',t(np)
c
cjk      print '("Highest frequency to be  computed: ",f7.0)', fmax
c      read *,fmax
cjk      print '("Resolution (i.e. points per peak): ",f7.0)', ofac
c      read *, ofac
c
      call AVEVAR(a,np,ave,var)
c
      print *,'                  average: ',ave
cjk      print *,' Variance of data points=',var
c
      hifac=fmax/fnyqav
      nf=0.5*ofac*hifac*np
cjk      print *, ' Frequency points to be output = ', nf

C SIZE THE FFT AS NEXT POWER OF 2 ABOVE NFREQT (MACC=4 interpolation points)
      nfreqt=ofac*hifac*np*4
      nfreq=64
      do while ( nfreq.LT.nfreqt )
        nfreq = nfreq*2
      enddo
      ndim=2*nfreq
cjk      print *, ' Workspace requirement (x3) =', ndim
      allocate (freq(ndim), STAT=ierror)
      if ( ierror /= 0 ) then
        write(*,*)"Error trying to allocate freq(ndim)"
        stop
      end if
      allocate (pow(ndim), STAT=ierror)
      if ( ierror /= 0 ) then
        write(*,*)"Error trying to allocate pow(ndim)"
        stop
      end if
      allocate (pha(ndim), STAT=ierror)
      if ( ierror /= 0 ) then
        write(*,*)"Error trying to allocate pha(ndim)"
        stop
      end if
c
      call POWER(t,a,np,ofac,hifac,freq,pow,pha,ndim,nfreq)
      do i=1,nfreq
        write (2,100) freq(i),sqrt(pow(i))
      enddo
100   format(1pe17.9,e14.6)
c
c Search Transform for peaks, print statistics for each peak found
c
c      dfreq= freq(2)-freq(1)
c      do i=2,nfreq-3
c        p1 = pow(i - 1)
c        p2 = pow(i)
c        p3 = pow(i + 1)
c        if ((p2.gt.p3).and.(p2.gt.p1)) then
c          alp2 = log(p1 / p2)
c          alp3 = log(p1 / p3)
c          bigp = alp2 / log(p2 / p3)
c          p = (1 + bigp)/(bigp - 1)/2.0
c          pp = p*p
c          f0 = freq(i) + p * dfreq
c          period = 1 / f0
c          p0 = p2 * exp(-alp3*p/4.0)
c          width = sqrt(dfreq**2 * bigp / (1 - bigp) / alp2)
c          phi1 = pha(i - 1)
c          phi2 = pha(i)
c          phi3 = pha(i + 1)
c          if ( phi1.lt.0 ) phi1 = phi1 + 2.0*pi
c          if ( phi2.lt.0 ) phi2 = phi2 + 2.0*pi
c          if ( phi3.lt.0 ) phi3 = phi3 + 2.0*pi
c          phi_0 = p*(p-1)*phi1/2 + (1-p*p)*phi2 + p*(p+1)*phi3/2;
c          if ( phi_0.gt.pi ) phi_0 = phi_0 - 2.0*pi
c          write (3,102) f0, sqrt(p0), period, phi_0
c102       format(1pe16.8,e12.4,e14.6,e14.6,e14.6)
c        endif
c      enddo
c
      write(*,*)
      write(*,*) '     FT save in file:     ft50.trf' 
      write(*,*)
      stop
      end
c
      subroutine power(x,y,n,ofac,hifac,freq,pow,pha,ndim,nout)
c
C GIVEN N DATA POINTS WITH ABSCISSAS X (WHICH NEED NOT BE EQUALLY
C SPACED) AND ORDINATES Y, AND GIVEN A DESIRED OVERSAMPLING FACTOR OFAC (A
C TYPICAL VALUE BEING 4 OR LARGER), THIS ROUTINE FILLS ARRAY freq WITH A
C SEQUENCE OF NOUT INCREASING FREQUENCIES (NOT ANGULAR FREQUENCIES) UP TO
C HIFAC TIMES THE "AVERAGE" NYQUIST FREQUENCY, AND FILLS ARRAY pow WITH
C THE VALUES OF THE POWER (deltaI/I squared) at those 
C FREQUENCIES.  THE ARRAYS X AND Y ARE NOT ALTERED. NWK, THE DIMENSION OF
C freq AND pow, MUST BE NDIM WHICH MUST BE A POWER OF 2
C
C adapted FROM PRESS AND RYBICKI, AP. J. 338, 277, as modified in 2nd Num Rec
C
      parameter (macc=4)
C MACC IS THE NUMBER OF INTERPOLATION POINTS PER 1/4 CYCLE OF HIGHEST FREQ.
      real x(n), y(n), freq(ndim), pow(ndim), pha(ndim)
      nout=0.5*ofac*hifac*n
C COMPUTE THE MEAN, VARIANCE, AND RANGE OF THE DATA
      CALL AVEVAR(Y,N,AVE,VAR)
      XMIN=X(1)
      XMAX=XMIN
      DO J=2,N
        IF (X(J).LT.XMIN)XMIN=X(J)
        IF (X(J).GT.XMAX)XMAX=X(J)
      ENDDO
      XDIF=XMAX-XMIN
C ZERO THE WORKSPACE
      DO J=1,NDIM
        freq(J)=0.
        pow(J)=0.
      ENDDO
      FAC=NDIM/(XDIF*OFAC)
      FNDIM=NDIM
C EXTIRPOLATE DATA INTO THE WORKSPACES
      DO J=1,N
        CK=1.+AMOD((X(J)-XMIN)*FAC,FNDIM)
        CKK=1.+AMOD(2.*(CK-1.),FNDIM)
        CALL SPREADD(Y(J)-AVE,freq,NDIM,CK,MACC)
        CALL SPREADD(1.,pow,NDIM,CKK,MACC)
      ENDDO
C TAKE THE FFT
      CALL REALFT(freq,NDIM,1)
      CALL REALFT(pow,NDIM,1)
      DF=1./(XDIF*OFAC)
c
      K=3
C compute the power for each frequency point
      DO J=1,NOUT
        HYPO=SQRT(pow(K)**2+pow(K+1)**2)
        pha(J)=ATAN2(pow(K+1),pow(K))
        HC2WT=0.5*pow(K)/HYPO
        HS2WT=0.5*pow(K+1)/HYPO
        CWT=SQRT(0.5+HC2WT)
        SWT=SIGN(SQRT(0.5-HC2WT),HS2WT)
        DEN=0.5*N+HC2WT*pow(K)+HS2WT*pow(K+1)
        CTERM=(CWT*freq(K)+SWT*freq(K+1))**2/DEN
        STERM=(CWT*freq(K+1)-SWT*freq(K))**2/(N-DEN)
        freq(J)=J*DF
c don't normalize power by twice the variance, but do normalize
c it by half the number of data points.
        pow(J)=(CTERM+STERM)/n*2.d0
        K=K+2
      ENDDO
      RETURN
      END
c
      SUBROUTINE REALFT(DATA,N,ISIGN)
C Calculates the Fourier transform of a set of N real-valued data
c points.  Replaces the data (which are stored in array DATA) by the
c positive frequency half of its complex Fourier transofrm.  The
c real-valued first and last components of the complex transform are
c returned as elements DATA(1) and DATA(2) respectively.  N must be a
c power of 2.  This routine also calculates the inverse transform of a
c complex data array if it is the transform of real data. (Result in this
c case must be multiplied by 2/N.)
C Numerical Recipes, p. 507.
      INTEGER ISIGN,N
      DIMENSION DATA(N)
      DOUBLE PRECISION THETA,WI,WPI,WPR,WR,WTEMP
      THETA=3.141592653589793d0/dble(n/2)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N/2,+1)
      ELSE
        C2=0.5
        THETA=-THETA
      ENDIF
      WPR=-2.0D0*SIN(0.5D0*THETA)**2
      WPI=SIN(THETA)
      WR=1.D0+WPR
      WI=WPI
      N2P3=N+3
      DO I=2,N/4
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
      ENDDO
      IF (ISIGN.EQ.1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
      ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N/2,-1)
      ENDIF
      RETURN
      END
c
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
c Replaces DATA by its discrete Fourier transform, if ISIGN is input as 1
C or replaces DATA by NN times its inverse discrete Fourier transifor
c transform if ISIGN=-1.  DATA is a complex array of length NN or, 
c equivalently, a real array of length 2*NN.  NN MUST be an integer
c power of 2 (this isn't checked for!)
      DOUBLE PRECISION WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(2*NN)
      N=2*NN
      J=1
      DO I=1,N,2
        IF (J.GT.I) THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
        DO WHILE ((M.GE.2).AND.(J.GT.M))
          J=J-M
          M=M/2
        ENDDO
        J=J+M
      ENDDO
      MMAX=2
      DO WHILE (N.GT.MMAX)
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO M=1,MMAX,2
          DO I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
          ENDDO
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
        ENDDO
        MMAX=ISTEP
      ENDDO
      RETURN
      END
C
      SUBROUTINE SPREADD(Y,YY,N,X,M)
C
C GIVEN ARRAY YY OF LENGTH N, EXTIRPOLATE (SPREAD) A VALUE Y INTO M
C ACTUAL ARRAY ELEMENTS THAT BEST APPROXIMATE THE "FICTIONAL" (I.E.
C POSSIBLY NON-INTEGER) ARRAY ELEMENT NUMBER X.  THE WEIGHTS USED ARE
C COEFFICIENTS OF THE LAGRANGE INTERPOLATING POLYNOMIAL.
C
      INTEGER m,n
      REAL x,y,yy(n)
      INTEGER ihi,ilo,ix,j,nden,nfac(10)
      REAL fac
      SAVE nfac
      DATA nfac /1,1,2,6,24,120,720,5040,40320,362880/
      IF (M.GT.10) PRINT *,'FACTORIAL TABLE TOO SMALL FOR SPREAD'
      IX=X
      IF (X.EQ.FLOAT(IX)) THEN
         YY(IX)=YY(IX)+Y
      ELSE
        ILO=MIN(MAX(INT(X-0.5*M+1.0),1),N-M+1)
        IHI=ILO+M-1
        NDEN=NFAC(M)
        FAC=X-ILO
        DO J=ILO+1,IHI
          FAC=FAC*(X-J)
        ENDDO
        YY(IHI)=YY(IHI)+Y*FAC/(NDEN*(X-IHI))
        DO J=IHI-1,ILO,-1
          NDEN=(NDEN/(J+1-ILO))*(J-IHI)
          YY(J)=YY(J)+Y*FAC/(NDEN*(X-J))
        ENDDO
      ENDIF
      RETURN
      END
c
      SUBROUTINE AVEVAR(DATA,N,AVE,VAR)
c Given array DATA of length N, returns its mean as AVE and its
c variance as VAR
c
      DIMENSION DATA(N)
      AVE=0.0
      VAR=0.0
      DO J=1,N
        AVE=AVE+DATA(J)
      ENDDO
      AVE=AVE/N
      DO J=1,N
        S=DATA(J)-AVE
        VAR=VAR+S*S
      ENDDO
      VAR=VAR/(N-1)
      RETURN
      END
