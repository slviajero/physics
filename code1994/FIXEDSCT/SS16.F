      PROGRAM SS16

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  ORGANISATION OF THE PROGRAM:                                        C
C    SUBROUTINE EDCS  CALCULATES ELASTIC DIFFERENTIAL CROSS SECTION    C
C    SUBROUTINE FORM  CALCULATES THE FORM FACTOR OF THE DISTRIBUTION   C
C    SUBROUTINE OXY   GENERATES A GAUSSIAN DENSITY DISTRIBUTION FOR    C
C                     OXYGEN-16                                        C
C    SUBROUTINE GAUSS GENERATES A GAUSSIAN DENSITY DISTRIBUTION FOR    C
C                     HELIUM-4                                         C
C    SUBROUTINE HOMO  GENERATES A HOMOHENEOUS DENSITY DISTRIBUTION     C
C    SUBROUTINE CONFI GENERATES CONFIGURATIONS AND CALCULATES THEIR    C
C                     WEIGHT;                                          C
C    REAL FUNCTION RAN CALCULATES A RANDOM NUMBER                      C
C    SUBROUTINE HCORE REJECTS CONFIGURATIONS WITH PARTICLES DISTANCES  C
C                     WITHIN A HARD CORE RADIUS                        C
C    SUBROUTINE OPTIC CALCULATES ELASTIC DIFFERENTIAL CROSS SECTION IN C
C                     OPTICAL MODEL APPROXIMATION                      C
C    COMPLEX FUNCTION VEFF  PROVIDES THE VALUES OF THE EFFECTIVE       C
C                     POTENTIAL FOR INTEGRATING THE SCHROEDINGER EQ.   C
C    SUBROUTINE BESSEL DETERMINES THE VALUES OF BESSEL AND NEUMANN     C
C                      FUNCTIONS                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PARAMETER USED IN THE WHOLE PROGRAM:                               C                          C
C      M      DIMENSIONALITY OF THE SPACE                              C
C      N      NUMBER OF PARTICLES                                      C
C      B      WIDTH OF THE HARMONIC OSCILLATOR WAVE FUNCTION           C
C   VARIABLES USED IN THE STOCHASTIC PART OF THE PROGRAM               C
C      NSTART  NUMBER OF INITIAL CONFIGURATIONS                        C
C      NKEEP   NUMBER OF CONFIGURATION KEPT AFTER EACH INITIAL         C
C              CONFIGURATION                                           C
C      X       INITIAL VALUE FOR RANDOM NUMBERS                        C
C      XD      COUNTER FOR RANDOM NUMBERS                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION X,XK

      COMMON/PARA/PI,NSTART,NKEEP
      COMMON/RG/X,XK

      PI=4.*ATAN(1.)

C  INITIALIZATION OF THE RANDOM NUMBER GENERATOR
      X=PI*1.D11
      XK=0.D00
      DO 10 K1=1,20
        TT=RAN()
10    CONTINUE
      NAUS=1

      OPEN (1,FILE='FORM',STATUS='UNKNOWN')
      OPEN (2,FILE='STOCHASTIC',STATUS='UNKNOWN')
      OPEN (3,FILE='CONFIG',STATUS='UNKNOWN')
      OPEN (4,FILE='OPTICAL',STATUS='UNKNOWN')
      OPEN (15,FILE='flucti',STATUS='UNKNOWN')
      OPEN (16,FILE='f1data',STATUS='UNKNOWN')
      open (17,FILE='forwar',status='UNKNOWN')

      REWIND(1)
      REWIND(2)
      REWIND(3)
      REWIND(4)
      REWIND(15)
      REWIND(16)
      rewind(17)



      READ(3,1010,ERR=100,END=100)NSTART,NKEEP
      PRINT*,' YOU ALREADY OWN A FILE WITH CONFIGURATIONS !   '
      PRINT*,'  '
      GOTO 200

100   JUMP=1
      WRITE(*,1000)
      WRITE(*,1001)
      WRITE(*,1002)
      READ(*,*,ERR=700)NANT
      IF (NANT.EQ.1) THEN
        CALL GAUSS
        NAUS=2
      ELSEIF (NANT.EQ.4) THEN
        CALL OPTIC
      ELSEIF (NANT.EQ.2) THEN
        CALL homo
        naus=2
      ELSEIF (NANT.EQ.3) THEN
        CALL OXY
        NAUS=2
      ELSE
        NAUS=0
      END IF
      IF (NAUS.EQ.0) GOTO 300
      IF (NAUS.NE.2) GOTO 100

200   JUMP=2
      REWIND(3)
      READ(3,1010)NSTART,NKEEP
      WRITE(*,1003)
      WRITE(*,1002)
      READ(*,*,ERR=700)NANT
      IF (NANT.EQ.1) THEN
        CALL FORM
      ELSEIF (NANT.EQ.2) THEN
        CALL EDCS_s
      ELSEIF (NANT.EQ.3) THEN
        CALL EDCS_sp
      ELSEIF (NANT.EQ.4) THEN
        CALL OPTIC
      ELSE
        NAUS=0
      END IF

      REWIND(3)
      READ(3,1010)NSTART,NKEEP
      IF (NAUS.NE.0) GOTO 200

300   CLOSE (1,STATUS='KEEP')
      CLOSE (2,STATUS='KEEP')
      CLOSE (3,STATUS='KEEP')
      CLOSE (4,STATUS='KEEP')
      CLOSE (15,STATUS='KEEP')
      CLOSE (16,STATUS='KEEP')
      CLOSE (17,STATUS='KEEP')

800   STOP

700   PRINT*,'GIVE THE LAST INPUT ONCE MORE!  '
      GOTO (100,200),JUMP

1000  FORMAT('  '/'  '/'   '/T12,' CHOOSE THE DESIRED ROUTINE OUT OF THE
     1 MENU'/'  '/'  '/T8,' 1 : GENERATION OF CONFIGURATIONS WITH GAUSSI
     2AN DISTRIBUTION'/t8,' 2 : generation of homogeneous distribution')
1001  FORMAT(T8,' 3 : GENERATION OF H.O. DISTRIBUTION FOR OXYGEN'/
     1T8,' 4 : OPTICAL POTENTIAL CALCULATION'
     2/T6,' ELSE  END THE PROGRAM'/'   ')
1003  FORMAT(T8,' 1 : CALCULATION OF THE FORM FACTOR'/T8,' 2 : CALCULATI
     1ON OF DIFFERENTIAL CROSS SECTION: s-wave'/T8,' 3 : calculation of
     2differential cross section: p-wave'/T8,' 4 : OPTICAL POTENTIAL CAL
     3CULATION'/T6,' ELSE  END THE PROGRAM'/'   ')
1002  FORMAT(T4,' TYPE THE NUMBER OF THE ROUTINE WHICH YOU WANT TO BE EX
     1ECUTED  ')
1010  FORMAT(2I8)

      END


      SUBROUTINE EDCS_s
      END


      SUBROUTINE EDCS_sp

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE CALCULATES THE SCATTERING AMPLITUDE AT 0 DEGREES IF    C
C  THE MOMENTUM OF THE INCOMING PARTICLE IS SMALLER THAN 0.2 INVERSE   C
C  FERMI; OTHERWISE DIFFERENTIAL CROSS SECTION. IN ADDITION, THE       C
C  NUMBER OF DIFFERENT CONFIGURATIONS, THE AVERAGE ACCEPTANCE          C
C  PROPABILITY, AND THE RMS-RADIUS IS DETERMINED.                      C                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PARAMETER:                                                         C
C      KAMAX   NUMBER OF VALUES OF THE SCATTERING ANGLE WHERE THE      C
C              DIFFERENTIAL CROSS SECTION IS CALCULATED                C
C      KMAX    MAXIMUM NUMBER OF P-WAVE AMPLITUDES READ FROM FILE      C
C              F1DATA                                                  C
C      Q0      UPPER LIMIT FOR CALCULATION OF SCATTERING AMPLITUDE IN  C
C              FORWARD DIRECTION ONLY                                  C
C      KF      NUMBER OF AVERAGED CROSS-SECTIONS FOR DETERMINING THE   C
C              MONTE CARLO ERRORS
C   IMPORTANT VARIABLES:                                               C
C      AM      CONTAINS THE MATRIX ELEMENTS OF THE MULTIPLE SCATTERING C
C              EQUATIONS                                               C
C      NTHET   ANGULAR STEP                                            C
C      QFACT   INCIDENT MOMENTUM STEP                                  C
C      FKKP    CONTAINS THE SCATTERING AMPLITUDES                      C
C      MULTI   MULTIPLICITY OF A CONFIGURATION                         C
C      PHI     INHOMOGENEOUS VECTOR OF THE MULTIPLE SCATTERING         C
C              EQUATIONS IN THE BEGINNING SOLUTION OF THESE EQUATIONS  C
C              IN THE END                                              C
C      PHI1    CONTAINS THE OUTGOING WAVE                              C
C      Q       MOMENTUM IN Z-DIRECTION                                 C
C      QANGLE  DIRECTIONS OF THE OUTGOING WAVE                         C
C      RMS     RMS-RADIUS OF THE DENSITY                               C
C      SLENG   PROJECTILE NUCLEON SCATTERING LENGTH                    C
C      SUM     TOTAL NUMBER OF CONFIGURATIONS                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  History (Stefan Lenz:
C  06.10.1991 : -comments added, loop variable names changed
C               -calculation of inhomogenous vector taken out
C                of matrix loops
C               -'if' statements taken out of matrix loop
C                 (replaced by do-loops)
C               -matrix routine prepared for sp-wave interaction
C               -s-wave part added to inhomogenous vector
C               -array dimensions set to sp-wave size
C               -NAG call changed to higher dimension
C
C  07.10.1991   -mechanism for multiple scattering lengthes
C                changed
C               -more loop variables renamed
C               -swave outgoing vector added
C               -output routines changed
C               -matrix elements added
C               -outgoing waves added
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(M=3,N=16,KAMAX=46,KMAX=10,Q0=.099999,kfm=10)
C
C lots of arrays (indices sometimes explained later)
C
      REAL FIELD(M,N),QANGLE(2,KAMAX),d(n,n,m+1)
     1     ,sum(kfm),rms(kfm)
C ------------------------------------ opt !
c    2     ,pangle(2,kamax),w1(kamax)
C-------------------------------------
      COMPLEX AM(4*N,4*N),PHI(4*N),FAKT,PSI1(4*N),PSI2(4*n)
     1  ,bm(4*n,4*n),ic,FKKP(kamax,kmax,kfm+1),PHI1(4*N,KAMAX),h(n,n,m)
     2  ,work(4*n)
C------------------------------------- opt!
c    2        ,phi2(4*n,kamax,kamax),y(kamax,kamax)
C-------------------------------------
      COMPLEX sleng(kmax),pleng(kmax)
C
C loop  variables
C
C   they count angles:
      INTEGER nangle,i
C   this one counts scattering lengthes
      INTEGER j
C   this one counts averaging steps
      INTEGER k
C   these count martices
      INTEGER n1,n2,n3

      COMMON/PARA/PI,NSTART,NKEEP
C
C  INPUT OF THE ANGULAR STEPS, MOMENTUM AND P-WAVE SCATTERING LENGTH
C  data is read from file 6 (f1data)
C
C   number of different scattering lengthes values
      READ(16,1045)KNUMB
C   array of s-wave scattering lengthes
        READ(16,*)SLENG
C   array of p-wave scattering lengthes
        READ(16,*)PLENG
C
      NTHET=180
      kangle=1
      mthet=1
C   read momentum and check range, for
C   very low momentum, only forward direction
C   is calculated
      READ(16,*)Q
      IF (Q.lt.Q0) GOTO 7
      READ(16,*)NTHET
      NTHET=MAX0(4,NTHET)
      mthet=nthet
      kangle=181
7     CONTINUE
C
C   kf is the number of averaging steps
C
      read(16,*)kf
      kf=min0(kf,kfm)
C
C  constants defined
C
C
      IC=CMPLX(0.,1.)
      sqr3=sqrt(3.)
      sqr45=sqrt(4.5)
C  number of configurations per averaging step:
      kmo=int(nstart*nkeep/(1.*kf))
C
C  CALCULATION OF UNIT VECTORS SPECIFYING THE DIRECTION OF THE
C  OUTGOING PROJECTILE
C  in qangle the first inex is the vector index (unlike d) !
C
      i=0
      DO 8 nangle=1,181,NTHET
        i=i+1
        QANGLE(1,i)=SIN(PI*(nangle-1)/180.)
        QANGLE(2,i)=COS(PI*(nangle-1)/180.)

C------------------------------------------------
c  zum test des optischen theorems
c       k11=0
c       do 11 k3=1,181,nthet
c         k11=k11+1
c         pangle(1,k11)=sin(pi*(k3-1)/90.)
c         pangle(2,k11)=cos(pi*(k3-1)/90.)
c         y(i,k11)=0.
c11       continue
C ------------------------------------------------

C
C scattering amplitudes set to zero
C FKKP(i,j,k) ist the scattering amplitude
C for the i th angle, j th scattering length
C and k th averaging step
C
        DO 8 K2=1,KNUMB
          do 8 k4=1,kf+1
            FKKP(i,K2,K4)=CMPLX(0.)
8     CONTINUE

      do 15 k=1,kf
        sum(k)=0.
        rms(k)=0.
15    continue
C
C central loop: NSTART*NKEEP configurations are loaded
C               FKKP is calculated and averaged over
C               kmo configurations (kf times)
C
      kfl=0
      DO 200 K1=1,NSTART*NKEEP
        if (mod(k1-1,kmo).eq.0) kfl=kfl+1
        if (kfl.gt.kf) kfl=kf
C load configuration:
        READ(3,*,END=800)FIELD
        READ(3,1010,END=800)RSQU,MULTI
C
C  BEGINNING OF THE CALCULATION OF THE SCATTERING AMPLITUDE FOR ONE
C  CONFIGURATION
C
        SUM(kfl)=SUM(kfl)+MULTI
        RMS(kfl)=RMS(kfl)+RSQU*MULTI
C
C  DETERMINATION OF THE PARTICLE DISTANCES
C  d(i,j,k) is the k th component of the distance
C  of i th and j th scattering center
C  d(i,j,4) is the length of this vector
C
         DO 20 K2=1,N-1
          do 10 k3=1,m
            d(k2,k2,k3)=0.
            h(k2,k2,k3)=cmplx(0.)
10        continue
          d(k2,k2,m+1)=0.
          DO 20 K3=K2+1,N
            R=0.
            DO 12 K4=1,M
              d(k2,k3,k4)=field(k4,k2)-field(k4,k3)
              R=R+d(k2,k3,k4)**2
12          CONTINUE
            r=sqrt(r)
            d(k2,k3,m+1)=r
            d(k3,k2,m+1)=r
            do 20 k4=1,m
              d(k2,k3,k4)=d(k2,k3,k4)/r
              d(k3,k2,k4)=-d(k2,k3,k4)
20      CONTINUE
        do 22 k2=1,m
          d(n,n,k2)=0.
          h(n,n,k2)=cmplx(0.)
22      continue
        d(n,n,m+1)=0.
C
C   CALCULATION OF HANKEL FUNCTIONS multiplied by q
C   h(i,j,k) is the hankelfunction of order k-1
C   with the distance between i and j as argument
C
        DO 25 k3=1,n-1
          do 25 k4=k3+1,n
            dd=d(k3,k4,m+1)
            h(k3,k4,1)=cexp(ic*q*dd)/dd
            h(k3,k4,2)=h(k3,k4,1)*cmplx(1./(q*dd),-1.)
            h(k3,k4,3)=h(k3,k4,2)*3./(q*dd)-h(k3,k4,1)
            h(k4,k3,1)=h(k3,k4,1)
            h(k4,k3,2)=h(k3,k4,2)
            h(k4,k3,3)=h(k3,k4,3)
25      continue
C
C   CALCULATION OF THE OUTGOING WAVES
C   this is simply exp(-i*k*x) (4times the same)
C
        DO 30 K3=1,4*n
          k5=mod(k3-1,n)+1
          K10=0
          DO 30 K4=1,KANGLE,MTHET
            K10=K10+1
            PHI1(K3,K10)=CEXP(-IC*Q*(FIELD(M,K5)*QANGLE(2,K10)+
     1                       FIELD(M-1,K5)*QANGLE(1,K10)))

C-----------------------------------------------------------
c  zum test des optischen theorems
c         k11=0
c         do 30 k6=1,kangle,mthet
c           k11=k11+1
c           PHI2(K3,K10,k11)=CEXP(-IC*Q*(FIELD(M,K5)*QANGLE(2,K10)+
c    1      FIELD(M-1,K5)*QANGLE(1,K10)*pangle(1,k11)+field(m-2,k5)
c    2      *qangle(1,k10)*pangle(2,k11)))
C-------------------------------------------------------------
30    CONTINUE
C
C   inhomogenous vector (the calculation is taken out
C   of the matrix loop).
C   as the incoming wave contains only components
C   in z-direction, only the Y(1,0) part is nonvanishing
C   this corresponds to the 1..N components of the vector
C
C   unlike in the pure p-wave program, the scattering length can not
C   betaken out so easily, thats why the final inhomogenous vector
C   is computed already here
C
C
C   the s-wave part is contained in the 3N+1..4N  components
C
      DO 31 n1=1,N
C      swave:
        PHI(n1+3*N)=cexp(ic*q*field(m,n1))
C      pwave:
        PHI(n1)    =sqr3*PHI(n1+3*N)
31    CONTINUE
      DO 32 n1=N+1,3*N
        PHI(n1)=0.
32    CONTINUE
C
C---------------------------------------------------------------------
C
C   CALCULATION OF THE MATRIX ELEMENTS: the first n elements correspond
c   to b(1,0,i), the next n to b(1,1,i), to b(1,-1,i)
C   the last to b(0,0,i)
C
C   am is the matrix which is to be inverted without the scattering
C   lengthes. the first index is a combination of (lmi) in
C   'boundary condition models' the second is a combination
C   of (lambda,mu,j)
C
C
c   calculation of matrix elements which couple b(1,0,i) to:
        DO 41 n1=1,N
C
c   b(1,0,j)
C
            do 411 n2=1,n
                if (n2.eq.n1) then
                  am(n1,n2)=1.
                else
                  am(n1,n2)=(-h(n1,n2,1)+h(n1,n2,3)*(3.*d(n1,n2,3)**2
     1                      -1.))
                end if
411         continue
C
c   b(1,1,j)
C
            do 412 n2=n+1,2*n
                if (n2-n.eq.n1) then
                  am(n1,n2)=0.
                else
                  am(n1,n2)=-sqr45*h(n1,n2-n,3)*d(n1,n2-n,3)*
     1               cmplx(d(n1,n2-n,1),d(n1,n2-n,2))
                end if
412         continue
C
c   b(1,-1,j)
C
            do 413 n2=2*n+1,3*n
                if (n2-2*n.eq.n1) then
                  am(n1,n2)=0.
                else
                  am(n1,n2)=sqr45*h(n1,n2-2*n,3)*d(n1,n2-2*n,3)*
     1                  cmplx(d(n1,n2-2*n,1),-d(n1,n2-2*n,2))
                end if
413         continue
C
C   b(0,0,j)
C
C
            do 414 n2=3*n+1,4*n
                if (n2-3*n.eq.n1) then
                  am(n1,n2)=0.
                else
                  am(n1,n2)=-ic*sqr3*d(n1,n2-3*n,3)*h(n1,n2-3*n,2)
                endif
414         continue
C
C end of b(1,0,i)
C
41      continue
C
C
c   calculation of matrix elements which couple b(1,1,i) to:
C
        do 42 n1=n+1,2*N
C
c   b(1,0,j)
C
            do 421 n2=1,n
                if (n2+n.eq.n1) then
                  am(n1,n2)=0.
                else
                  am(n1,n2)=-sqr45*h(n1-n,n2,3)*d(n1-n,n2,3)*
     1                    cmplx(d(n1-n,n2,1),-d(n1-n,n2,2))
                end if
421         continue
C
c   b(1,1,j)
C
            do 422 n2=n+1,2*n
                if (n2.eq.n1) then
                  am(n1,n2)=1.
                else
                  am(n1,n2)=(-h(n1-n,n2-n,1)-h(n1-n,n2-n,3)*(3*
     1                    d(n1-n,n2-n,3)**2-1.)*.5)
              end if
422         continue
C
c   b(1,-1,j)
C
            do 423 n2=2*n+1,3*n
                if (n2-n.eq.n1) then
                  am(n1,n2)=0.
                else
                  am(n1,n2)=-1.5*h(n1-n,n2-2*n,3)*
     1                 cmplx(d(n1-n,n2-2*n,1)**2-d(n1-n,n2-2*n,2)**2,
     2                 -2.*d(n1-n,n2-2*n,1)*d(n1-n,n2-2*n,2))
                end if
423         continue
C
C   b(0,0,j)
C
C
            do 424 n2=3*n+1,4*n
                if(n2-2*n.eq.n1) then
                   am(n1,n2)=0.
                else
                   am(n1,n2)=ic*sqr45/sqr3*
     1                cmplx(d(n1-n,n2-3*n,1),-d(n1-n,n2-3*n,2))*
     2                h(n1-n,n2-3*n,2)
                endif
424         continue
C
C   end of b(1,1,i)
C
42      continue
C
c   calculation of matrix elements which couple b(1,-1,i) to:
C
          do 43 n1=2*n+1,3*n
C
c   b(1,0,j)
C
            do 431 n2=1,n
                if (n2+2*n.eq.n1) then
                  am(n1,n2)=0.
                else
                  am(n1,n2)=sqr45*h(n1-2*n,n2,3)*d(n1-2*n,n2,3)*
     1                 cmplx(d(n1-2*n,n2,1),d(n1-2*n,n2,2))
                end if
431         continue
C
c   b(1,1,j)
C
            do 432 n2=n+1,2*n
                if (n2+n.eq.n1) then
                  am(n1,n2)=0.
                else
                  am(n1,n2)=-1.5*h(n1-2*n,n2-n,3)*
     1                cmplx(d(n1-2*n,n2-n,1)**2-d(n1-2*n,n2-n,2)**2,
     2                2.*d(n1-2*n,n2-n,1)*d(n1-2*n,n2-n,2))
                end if
432         continue
C
c   b(1,-1,j)
C
            do 433 n2=2*n+1,3*n
                if (n2.eq.n1) then
                    am(n1,n2)=1.
                else
                  am(n1,n2)=(-h(n1-2*n,n2-2*n,1)-h(n1-2*n,n2-2*n,3)*
     1                 (3.*d(n1-2*n,n2-2*n,3)**2-1.)*.5)
                end if
433         continue
C
C   b(0,0,j)
C
C
            do 434 n2=3*n+1,4*n
                if(n2-n.eq.n1) then
                   am(n1,n2)=0.
                else
                   am(n1,n2)=-ic*sqr45/sqr3*
     1                cmplx(d(n1-2*n,n2-3*n,1),+d(n1-2*n,n2-3*n,2))*
     2                h(n1-2*n,n2-3*n,2)
                endif
434         continue
C
C  end of b(1,-1,i)
C
43        continue
C
c   calculation of matrix elements which couple b(0,0,i) to:
C
          do 44 n1=3*n+1,4*n
C
c   b(1,0,j)
C
            do 441 n2=1,n
                if (n2+3*n.eq.n1) then
                  am(n1,n2)=0.
                else
                   am(n1,n2)=-ic*sqr3*
     1                d(n1-3*n,n2,3)*
     2                h(n1-3*n,n2,2)
                end if
441         continue
C
c   b(1,1,j)
C
            do 442 n2=n+1,2*n
                if (n2+2*n.eq.n1) then
                  am(n1,n2)=0.
                else
                   am(n1,n2)=ic*sqr45/sqr3*
     1                cmplx(d(n1-3*n,n2-n,1),+d(n1-3*n,n2-n,2))*
     2                h(n1-3*n,n2-n,2)
                end if
442         continue
C
c   b(1,-1,j)
C
            do 443 n2=2*n+1,3*n
                if (n2+n.eq.n1) then
                   am(n1,n2)=0.
                else
                   am(n1,n2)=-ic*sqr45/sqr3*
     1                cmplx(d(n1-3*n,n2-2*n,1),-d(n1-3*n,n2-2*n,2))*
     2                h(n1-3*n,n2-2*n,2)
                end if
443         continue
C
C   b(0,0,j)
C
C
            do 444 n2=3*n+1,4*n
                if (n2.eq.n1) then
                   am(n1,n2)=1.
                else
                   am(n1,n2)=-h(n1-3*n,n2-3*n,1)
                endif
444         continue
C
C  end of b(0,0,i)
C
44        continue
C
C  end of matrix build-up
C----------------------------------------------------------


C
C  LOOP FOR CALCULATING THE SCATTERING AMPLITUDE FOR EACH CONFIGURATION
C  AND ALL DIFFERENT P-WAVE AMPLITUDES
C
C  here things become complicated, because the scattering length can
C  to be taken out like in the pure p-wave case
C  thus the matrix am has to be multiplied with the appropriate
C  factors when transscribed to bm
C
      DO 200 j=1,KNUMB
C
C-------------------------------------------------------------
C  here the matrix am is transscribed to bm and multiplied
C  by the scattering lengthes out of PLENG and SLENG
C  this is not so efficient as in the p-wave case, but
C  the advantage of computing the sperical harmonics
C  once for many scattering lengthes remains
C
C
C the first 3*N rows are multiplied with PLENG
C from 3*N+1 to 4*N : SLENG
C keep in mind, that n2 is associated with lambda!
        DO 51 n1=1,4*N
          do 50 n2=1,3*N
            if (n1.ne.n2) then
               bm(n1,n2)=am(n1,n2)*PLENG(j)
            else
               bm(n1,n2)=am(n1,n2)
            endif
50        CONTINUE
          do 51 n2=3*N+1,4*N
            if (n1.ne.n2) then
               bm(n1,n2)=am(n1,n2)*SLENG(j)
            else
               bm(n1,n2)=am(n1,n2)
            endif
51      CONTINUE
C--------------------------------------------------------------
C
C preparation of the NAG-lib call
C
C       write(1,*)'---------ss16---------------'
C       write(1,*)'Confi'
C       write(1,*)field
C       write(1,*)'pre-matrix'
C       write(1,*)am
C       write(1,*)'matrix'
C       write(1,*)bm
C       write(1,*)'inhomo'

C       do 7777 iii=1,n
C         do 7777 jjj=1,n
C7777       write(1,*)bm(iii+3*n,jjj+3*n)



        ifail=0
        m0=1
        m1=4*n
        call f04adf(bm,m1,phi,m1,m1,m0,psi1,m1,work,ifail)

C       write(1,*)'loesung'
C       write(1,*)psi1



C
C the result of the calculation is now in psi1
C
C---------------------------------------------------------------
C
C   AVERAGING OF THE INDIVIDUAL SCATTERING AMPLITUDES
C   here the result vector psi1 is projected on
C   the outgoing waves
C
          DO 200 n1=1,4*N
            i=0
            DO 200 nangle=1,KANGLE,MTHET
              i=i+1
C-----------------------------------
c  zum test des optischen theorems
c         k11=0
c         do 200 k9=1,kangle,mthet
c           k11=k11+1
C-----------------------------------
            if (n1.le.n) then
C-----------------------------------
c              if (k11.eq.1) then
C-----------------------------------
C
C  Y(1,0) part of outgoing wave:
C
              FKKP(i,j,Kfl)=FKKP(i,j,Kfl)+MULTI*PsI1(n1)*
     1             PHI1(n1,i)*qangle(2,i)*sqr3*PLENG(j)
C-----------------------------------
c             end if
c             y(i,k11)=psi1(n1)*phi2(n1,i,k11)*qangle(2,i)*sqr3
c    1                    *PLENG(j)+y(i,k11)
C-----------------------------------
            else if ((n1.gt.n).and.(n1.le.2*n)) then
C-----------------------------------
c             if (k11.eq.1) then
C-----------------------------------
C
C  Y(1,1) part of outgoing wave:
C
              FKKP(i,j,Kfl)=FKKP(i,j,Kfl)+MULTI*PsI1(n1)*
     1        PHI1(n1,i)*cmplx(0.,-qangle(1,i)*sqr45/sqr3)*PLENG(j)
C------------------------------------
c             end if
c             y(i,k11)=psi1(n1)*phi2(n1,i,k11)*qangle(1,i)*
c    1        cmplx(-pangle(2,k11),-pangle(1,k11))*sqr45/sqr3*PLENG(k)
C    2                       +y(i,k11)
C------------------------------------
            else if ((n1.gt.2*n).and.(n1.le.3*n)) then
C------------------------------------
c             if (k11.eq.1) then
C------------------------------------
C
C  Y(1,-1) part of outgoing wave:
C
              FKKP(i,j,Kfl)=FKKP(i,j,Kfl)+MULTI*PsI1(n1)*
     1        PHI1(n1,i)*cmplx(0.,-qangle(1,i)*sqr45/sqr3)*PLENG(j)
C------------------------------------
c             end if
c             y(i,k11)=psi1(n1)*phi2(n1,i,k11)*qangle(1,i)*
c    1        cmplx(pangle(2,k11),-pangle(1,k11))*sqr45/sqr3*PLENG(j)
C    2                +y(i,k11)
C------------------------------------
            else if (n1.gt.3*n) then
C
C  Y(0,0) part of outgoing wave:
C
              FKKP(i,j,Kfl)=FKKP(i,j,Kfl)+MULTI*PsI1(n1)*
     1        PHI1(n1,i)*SLENG(j)
            end if
200   CONTINUE
C
C  after termination of this loop, FKKP(i,j,kfl) contains
C  the averaged scattering amplitude of
C     angle i
C     scattering length j
C     averaging step kfl
C
C---------------------------------------------------------------
C
C  STORING OF SCATTERING ANGLES, AVERAGED SCATTERING AMPLITUDES AND
C  DIFFERENTIAL CROSS SECTION
C
C  here the information is stored to various files
C     FILE 2 : STOCHSATIC  (general result file)
C     FILE 5 : flucti      (data in PLOTable form)
C     FILE 7 : forwar      (contains forward scattering information)
C

      DO 300 K1=1,KNUMB
        k10=0
        do 220 k2=1,kangle,mthet
          k10=k10+1
          qangle(1,k10)=0.
          qangle(2,k10)=0.
220     continue

C
C  start message and scattering length written to STOCHASTIC
C
        WRITE(2,1020)SLENG(K1),PLENG(K1)

        summe=0.
        rrms=0.
C
C  start message written to flucti
C
        write(15,1015)nstart*nkeep,kf,q,sleng(k1),pleng(K1)
        do 210 k2=1,kf
          summe=summe+sum(k2)
          rrms=rrms+rms(k2)
          k10=0
          kk=k2*kmo
          if (k2.eq.kf) kk=nstart*nkeep
          if (mod(k2,2).eq.0) write(15,1025)1,kk,summe,(rrms/(summe*
     1         n)+.64)**.5
          do 210 k3=1,kangle,mthet
            k10=k10+1
            fkkp(k10,k1,kf+1)=fkkp(k10,k1,kf+1)+fkkp(k10,k1,k2)
            if (mod(k2,2).eq.0) then
              write(15,1000)k2,k3-1.,10.*cabs(fkkp(k10,k1,kf+1)
     1             /summe)**2.
            end if
            hilf=10.*cabs(fkkp(k10,k1,k2)/sum(k2))**2./kf
            qangle(1,k10)=qangle(1,k10)+hilf**2.
            qangle(2,k10)=qangle(2,k10)+hilf
210     continue

        k10=0
C
C here the forwrd scattering amplitude is writen to
C 'forwar' as a function of the scattering length
C
        write(17,1000)98,cabs(pleng(k1)),real(fkkp(1,k1,kf+1))/summe
        write(17,1000)99,cabs(pleng(k1)),aimag(fkkp(1,k1,kf+1))/summe
        do 230 k2=1,kangle,mthet
          k10=k10+1
          write(15,1035)kf+1,k2-1.,10.*cabs(fkkp(k10,k1,kf+1)/
     1    summe)**2.,sqrt(qangle(1,k10)-qangle(2,k10)**2./kf),0.
230     continue

        WRITE(2,1030)Q
        IF (K1.EQ.1) WRITE(2,1050)
        K10=0
        DO 250 k2=1,KANGLE,MTHET
          K10=K10+1

c  zum test des optischen theorems
c         w1(k10)=0.
            WRITE(2,1040)k2-1,FKKP(k10,K1,Kf+1)/SUMme,10.*
     1            (CABS(FKKP(k10,K1,Kf+1))/SUMme)**2

c           k11=-1
c           do 240 k9=1,kangle-2*mthet,2*mthet
c             k11=k11+2
c             w1(k10)=w1(k10)+(cabs(y(k10,k11))**2+4.*cabs(
c    1             y(k10,k11+1))**2+cabs(y(k10,k11+2)*)**2)/3.
c240         continue
c            w1(k10)=w1(k10)*qangle(1,k10)
250       continue
c          zz=0.
c         k10=-1
c         do 260 k9=1,kangle-2*mthet,2*mthet
c           k10=k10+2
c           zz=zz+(w1(k10)+4.*w1(k10+1)+w1(k10+2))*2.*(pi*
c    1           mthet)**2/(3.*(kangle-1.)**2)
c260     continue
c       write(2,1070)q,4.*pi*aimag(y(1,1))/q,zz
300   continue

C   STORING OF INFORMATION CONCERNING THE CONFIGURATIONS AND THE RMS; THE
C   .64 BEING ADDED IS THE SQUARE OF THE NUCLEON RADIUS
      WRITE(2,1060)NSTART*NKEEP,NSTART*NKEEP/SUMme,(RrMS/(SUMme*N)+
     1   .64)**.5

800   RETURN

1000  format(i3,2e12.4)
1015  format(t5,'total # of different konfig. ',i8,t45,'fluktuation para
     1meter ',i3/t5,'momentum ',e12.4,' s-wave amplitude ',2e12.4,
     2' p-wave amplitude ',2e12.4)
1025  format(t5,'interval',i8,'--',i8/' number of konfigurations ',e12.4
     1,' rms-radius ',e12.4)
1035  format(i3,4e12.4)
1005  FORMAT(E12.4)
1010  FORMAT(E12.4,I8)
1020  FORMAT(T5,'STOCHASTIC CALCULATION:'/'    '
     1/T3,'SCATTERING LENGTH [FM]: s:',E12.4,' + I* ',E12.4,'  p:',
     2E12.4,' + I* ',E12.4,'  ')
1030  FORMAT(T3,'MOMENTUM [FM**(-1)]:   ',E12.4/'  ')
1045  FORMAT(I3)
1050  FORMAT(T22,'SCATTERING AMPLITUDE [FM]'/T3,
     1'SCATTERING ANGLE',T22,'REAL PART',T33,'IMAGINARY PART',T52,'CROSS
     2 SECTION [mb]')
1070  format(t3,'test of optical theorem'/t3,'momentum',t25,'Im F(0)',
     1t50,'F(Omega)dOmega'/e12.4,t22,e12.4,t50,e12.4)
1060  FORMAT(' '/T5,'NUMBER OF DIFFERENT CONFIGURATIONS',T60,I8/T5,
     1'AVERAGE ACCEPTANCE PROPABILITY',T60,E12.4/T5,'RMS-RADIUS',T60,
     2E12.4/'  ')
1040  FORMAT(T8,I4,T20,2E12.4,T52,E12.4)
1090  FORMAT(T5,' # OF VALUES FOR THE INCIDENT MOMENTUM K ( <= ',I2,')
     1 ')
1100  FORMAT(I3,2E12.4)
1110  FORMAT(I3,4E12.4)

      END






      SUBROUTINE FORM

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS SUBROUTINE CALCULATES THE FORM FACTOR AND ITS FLUCTUATION     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PARAMETER:                                                         C
C      MAXI  MAXIMUM NUMBER OF VALUES OF MOMENTUM TRANSFER             C
C   IMPORTANT VARIABLES:                                               C
C      FACTOR  MOMENTUM TRANSFER STEP                                  C
C      FLUCT   CONTAINS THE FLUCTUATIONS                               C
C    L/UBOUND  LOWER/UPPER BOUND OF MOMENTUM TRANSFER                  C
C      MULTI   MULTIPLICITY OF A CONFIGURATION                         C
C      Q       MOMENTUM TRANSFER                                       C
C      REFORM  FORM FACTOR FOR INCOMING WAVES EVALUATED FOR MOMENTUM   C
C              TRANSFER ALONG THE 3 COORDINATE AXES                    C
C      RMS     RMS-RADIUS OF THE DENSITY                               C
C      SUM     TOTAL NUMBER OF CONFIGURATIONS                          C
C     W/W1     AUXILIARY FIELDS FOR EVALUATION OF FLUCTUATIONS         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(M=3,N=16,B=1.7,MAXI=21)
      REAL W(MAXI),W1(MAXI),LBOUND,REFORM(M,MAXI),FIELD(M,N),FLUCT(MAXI)

      COMMON/PARA/PI,NSTART,NKEEP

1     JUMP=1
      WRITE(*,1090)MAXI
      READ(*,*,ERR=700)MOMENT
      MOMENT=MIN0(MOMENT,MAXI)
      IF (MOMENT.LT.1) MOMENT=1
      JUMP=2
2     PRINT*,' LOWER LIMIT OF Q [FM**(-1)]?  '
      READ(*,*,ERR=700)LBOUND
      JUMP=3
3     PRINT*,' UPPER LIMIT OF Q [FM**(-1)]?  '
      READ(*,*,ERR=700)UBOUND

      DO 5 K1=1,MOMENT
        FLUCT(K1)=0.
        W(K1)=0.
        DO 5 K2=1,M
          REFORM(K2,K1)=0.
5     CONTINUE
      RMS=0.
      SUM=0.

      IF (MOMENT.LE.1) THEN
        FACTOR=0.
      ELSE
        FACTOR=(UBOUND-LBOUND)/(MOMENT-1)
      END IF

      DO 40 J3=1,NSTART
        DO 25 J0=1,MOMENT
          W1(J0)=0.
25      CONTINUE
        SUM1=0.
        DO 30 J4=1,NKEEP
          READ(3,*,END=800)FIELD
          READ(3,1015,END=800)RSQU,MULTI

C  BEGINNING OF THE CALCULATION OF THE FORMFACTOR FOR ONE CONFIGURATION
          RMS=RMS+MULTI*RSQU
          DO 20 NQ=1,MOMENT
            Q=LBOUND+(NQ-1)*FACTOR
            DO 20 J2=1,N
              W1(NQ)=W1(NQ)+MULTI*COS(Q*FIELD(3,J2))/N
              DO 20 J1=1,M
                REFORM(J1,NQ)=REFORM(J1,NQ)+MULTI*COS(Q*FIELD(J1,J2))/N
20        CONTINUE
          SUM=SUM+MULTI
          SUM1=SUM1+MULTI
30      CONTINUE
        DO 40 NQ=1,MOMENT
          FLUCT(NQ)=FLUCT(NQ)+(W1(NQ)/SUM1)**2
          W(NQ)=W(NQ)+W1(NQ)/SUM1
40    CONTINUE

C  STORING VALUES OF MOMENTUM TRANSFER, FORMFACTORS AND FLUCTUATIONS
      WRITE(1,1010)
      DO 50 NQ=1,MOMENT
        FLUCT(NQ)=FLUCT(NQ)/NSTART-(W(NQ)/NSTART)**2
        Q=LBOUND+(NQ-1)*FACTOR
        WRITE(1,1020)Q,(REFORM(I1,NQ)/SUM,I1=1,M),(FLUCT(NQ)/NSTART)**.5
50    CONTINUE

C   STORING OF INFORMATION CONCERNING THE CONFIGURATIONS AND THE RMS; THE
C   .64 BEING ADDED IS THE SQUARE OF THE NUCLEON RADIUS
      WRITE(1,1030)NSTART*NKEEP,NSTART*NKEEP/SUM,(.64+RMS/(SUM*N))**.5

800   RETURN

700   PRINT*,'GIVE THE LAST INPUT ONCE MORE!  '
      GOTO (1,2,3),JUMP

1005  FORMAT(E12.4)
1015  FORMAT(E12.4,I8)
1010  FORMAT(' FORMFACTOR:'/T22,'  MOMENTUM TRANSFER ALONG:'/T8,'MOMENTU
     1M',T22,'X-DIRECTION',T34,'Y-DIRECTION',T46,'Z-DIRECTION',T64,'FLUC
     2TUATION')
1020  FORMAT(T5,E12.4,T20,3E12.4,T62,E12.4)
1030  FORMAT('  '/T3,' NUMBER OF DIFFERENT CONFIGURATIONS: ',T60,I8/T3,
     1' AVERAGE ACCEPTANCE PROPABILITY',T60,E12.4/T4,'RMS-RADIUS',T60,
     2E12.4)
1110  FORMAT(I3,2E12.4)
1120  FORMAT(I3,4E12.4)
1090  FORMAT(T5,' # OF VALUES OF MOMENTUM TRANSFER Q FOR WHICH THE FOURI
     1ER TRANSFORM'/T5,' IS CALCULATED (MAXIMUM ',I2,') ')

      END


      SUBROUTINE OXY

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS ROUTINE GENERATES CONFIGURATIONS WITH GAUSSIAN DISTRIBUTION   C
C   WITHOUT CENTER OF MASS CORRELATION OF THE NUCLEONS BY DIRECT       C
C   SAMPLING OF ALL POSITIONS                                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(M=3,N=16,B=1.7)

      real field(M,N)

      COMMON/PARA/PI,NSTART,NKEEP
      COMMON/HARDC/CORE

      print*,'WIEVIELE KONFIGURATIONEN ?  '
      read*,NSTART
      WRITE(3,1000)NSTART,1
      PRINT*,'HARD-CORE RADIUS [FM] ?                       '
      READ(*,*)CORE
      IF (CORE.GT.0.) CORE=CORE*CORE
      PRINT*,' WHICH FACTOR OF EXTENSION ? '
      READ*,FF
      NKEEP=1
      r0=b*FF

      do 500 k1=1,NSTART

C  GENERATION OF 4 PARTICLES IN S-WAVE-WAVEFUNCTION
15      do 20 k2=1,2
          do 20 k3=1,3
            phi=2*pi*ran()
            r=b*sqrt(-alog(ran()))
            field(k3,k2)=r*cos(phi)
            field(k3,k2+2)=r*sin(phi)
20      continue

C  GENERATION OF 8 PARTICLES IN PX- AND PY- WAVEFUNCTIONS
        do 30 k2=1,4
          phi=2*pi*ran()
          r=b*sqrt(-alog(ran()))
          field(3,4+k2)=r*cos(phi)
          field(3,8+k2)=r*sin(phi)
30      continue
        do 40 k2=1,8
          phi=2*pi*ran()
35        r=b*sqrt(-alog(ran()))
            rp=r0*SQRT(ran())
          if (rp.gt.r) goto 35
          field(1,4+k2)=r*cos(phi)
          field(2,4+k2)=r*sin(phi)
40      continue

C  GENERATION OF 4 PARTICLES IN PZ-WAVEFUNCTION
        do 50 k2=1,4
          phi=2.*pi*ran()
          r=b*sqrt(-alog(ran()))
          field(1,12+k2)=r*cos(phi)
          field(2,12+k2)=r*sin(phi)
45        r=b*sqrt(-alog(ran()))
            rp=r0*ran()
          if (rp.gt.r) goto 45
          field(3,12+k2)=r
50      continue
        CALL HCORE(*15,FIELD)
        write(3,*)field

c  CALCULATION OF THE SQUARE OF THE EXTENSION OF THE PARTICLES
C  AND STORATION
        r=0.
        do 60 k2=1,16
          do 60 k3=1,3
            r=r+field(k3,k2)**2
60      continue
        write(3,1010)r,1
500   continue

1000  format(2i8)
1010  FORMAT(E12.4,I8)
      RETURN
      end


      SUBROUTINE GAUSS

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS ROUTINE GENERATES CONFIGURATIONS WITH GAUSSIAN DISTRIBUTION   C
C   WITH CENTER OF MASS CORRELATION OF THE NUCLEONS BY USE OF THE      C
C   METROPOLIS ALGORITHM                                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   IMPORTANT VARIABLES:                                               C
C      AFIELD CONTAINES ACCEPTED CONFIGURATION                         C
C      C      WIDTH PARAMETER OF AUXILIARY DENSITY "SIGMA"             C
C      CORE   HARD CORE RADIUS                                         C
C      FIELD  CONTAINS TRIAL CONFIGURATION                             C
C      MULTI  MULTIPLICITY OF THE CONFIGURATION                        C
C      NCOUNT COUNTS FOR ACCEPTABLE CONFIGURATIONS                     C
C      NKEEP1 COUNTS THE CONFIGURATIONS STORED FOR AN INITIAL          C
C             CONFIGURATION                                            C
C      RAND   CONTAINS THE RANDOM NUMBERS CORRESPONDING TO THE         C
C             ACCEPTED AND TRIAL CONFIGURATION                         C
C      R      WEIGHT OF THE ACCEPTED AND TRIAL CONFIGURATION           C
C      RSQU   SUM OVER SQUARES OF PARTICLE DISTANCES                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(M=3,N=16,B=1.7)

      REAL RSQU(2),RAND(M,N-1,2),FIELD(M,N),AFIELD(M,N),R(2)

      COMMON/PARA/PI,NSTART,NKEEP/OLD/FIELD,RAND,R,RSQU
      COMMON/PARA1/C,WIDTH/HARDC/CORE

1     JUMP=1
      REWIND(3)
      PRINT*,'# OF INITIAL CONFIGURATIONS                              '
      READ(*,*,ERR=700)NSTART
      NSTART=IABS(NSTART)
      JUMP=2
2     PRINT*,'# OF CONFIGURATIONS GENERATED BEFORE ACCEPTING THE FIRST '
      READ(*,*,ERR=700)NNOT
      JUMP=3
3     PRINT*,'# OF CONFIGURATIONS KEPT FOR EACH INITIAL CONFIGURATION  '
      READ(*,*,ERR=700)NKEEP
      NKEEP=IABS(NKEEP)
      JUMP=4
4     PRINT*,'# OF CONFIGURATIONS LEFT OUT BETWEEN TWO ACCEPTED CONFIGUR
     1ATIONS '
      READ(*,*,ERR=700)NJUMP
      JUMP=5
5     PRINT*,'JUMP WIDTH? RECOMMENDED VALUE IS: WIDTH=0.2   '
      READ(*,*,ERR=700)WIDTH
      WIDTH=ABS(WIDTH)
      WIDTH=AMIN1(WIDTH,.5)
      JUMP=6
6     PRINT*,'HARD-CORE RADIUS [FM] ?                       '
      READ(*,*,ERR=700)CORE
      IF (CORE.GT.0.) CORE=CORE*CORE

      C=B*SQRT(3.)/2.

C  STORING THE NUMBER OF INITIAL CONFIGURATIONS AND NUMBER OF
C  CONFIGURATIONS GENERATED FROM INITIAL ONE
      WRITE(3,1000)NSTART,NKEEP
      DO 50 K1=1,NSTART

C   CALCULATION OF INITIAL CONFIGURATION
        CALL CONFI(1)
        DO 8 K2=1,M
          DO 8 K3=1,N
            AFIELD(K2,K3)=FIELD(K2,K3)
8       CONTINUE
        NKEEP1=0
        NCOUNT=0
10      NCOUNT=NCOUNT+1
          MULTI=0
20        MULTI=MULTI+1

C   UPDATING OF CONFIGURATIONS
          CALL CONFI(2)
          TT=RAN()
          IF (R(2)/R(1).LT.RAN()) GOTO 20

          IF (((NKEEP1.EQ.0).AND.(NCOUNT.LT.NNOT)).OR.(NCOUNT.LT.NJUMP))
     1       GOTO 30

C  SELECTION OF CONFIGURATIONS WITH PARTICLES CLOSER THAN THE HARD CORE
C  RADIUS
          IF (CORE.GT.0.) CALL HCORE(*30,AFIELD)

C  STORING CONFIGURATION, PARAMETERS AND READJUSTMENT OF COUNTERS
          WRITE(3,*)AFIELD
          WRITE(3,1005)RSQU(1),MULTI
          NCOUNT=0
          NKEEP1=NKEEP1+1

C  REPLACEMENT OF THE OLD BY THE NEW CONFIGURATION
30        DO 40 K4=1,M
            AFIELD(K4,N)=FIELD(K4,N)
            DO 40 K3=1,N-1
              AFIELD(K4,K3)=FIELD(K4,K3)
              RAND(K4,K3,1)=RAND(K4,K3,2)
40        CONTINUE
          R(1)=R(2)
          RSQU(1)=RSQU(2)
        IF (NKEEP1.LT.NKEEP) GOTO 10
50    CONTINUE

800   RETURN

700   PRINT*,'GIVE THE LAST INPUT ONCE MORE!  '
      GOTO (1,2,3,4,5),JUMP

1000  FORMAT(2I8)
1005  FORMAT(E12.4,I8)
1010  FORMAT(E12.4)
      END



      SUBROUTINE CONFI(I)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   AN INITIAL CONFIGURATION IS GENERATED FOR I=1, A TRIAL ONE FOR I=2 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   IMPORTANT VARIABLES:                                               C
C      R      WEIGHT OF THE ACCEPTED/TRIAL CONFIGURATION               C
C      RAND   ARRAY OF RANDOM NUMBERS SPECIFIING THE POSITIONS         C
C      Y      EXPONENT APPEARING IN THE WHEIGHT FUNCTION               C
C      R0     RADIAL DISTANCE FROM THE CENTER                          C
C      COTHET COSINE OF THE POLAR ANGLE                                C
C      SITHET SINE OF THE POLAR ANGLE                                  C
C      PHI    AZIMUTHAL ANGLE                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(M=3,N=16,B=1.7)

      REAL RAND(M,N-1,2),FIELD(M,N),R(2),RSQU(2)

      COMMON/PARA/PI,NSTART,NKEEP/OLD/FIELD,RAND,R,RSQU
      COMMON/PARA1/C,WIDTH/HARDC/CORE

5     R(I)=1.
      Y=0.
      RSQU(I)=0.

C  INITIALISATION OF A SET OF RANDOM NUMBERS AS STARTING VECTORS FOR
C  THE RANDOM WALK
      DO 10 K1=1,M
        FIELD(K1,N)=0.
        DO 10 K2=1,N-1
          IF (I.EQ.1) THEN
            RAND(K1,K2,I)=RAN()
          ELSE
            RAND(K1,K2,2)=RAND(K1,K2,1)+2.*WIDTH*(RAN()-.5)
            IF (RAND(K1,K2,2).GT.1.) RAND(K1,K2,2)=RAND(K1,K2,2)-1.
            IF (RAND(K1,K2,2).LT.0.) RAND(K1,K2,2)=RAND(K1,K2,2)+1.
          END IF
10    CONTINUE

C  CALCULATION OF THE PARTICLE POSITIONS
      DO 20 K1=1,N-1
        R0=C*SQRT(-ALOG(RAND(1,K1,I)))
        COTHET=2.*RAND(2,K1,I)-1.
        SITHET=SQRT(1.-COTHET*COTHET)
        PHI=2.*PI*RAND(3,K1,I)
        R(I)=R(I)*R0
        RSQU(I)=RSQU(I)+R0*R0
        Y=Y-R0*R0/3.
        FIELD(1,K1)=R0*SITHET*COS(PHI)
        FIELD(2,K1)=R0*SITHET*SIN(PHI)
        FIELD(3,K1)=R0*COTHET

C   LAST POSITION DETERMINED BY THE DEMAND FOR A FIXED CENTER OF MASS
        DO 20 K2=1,M
          FIELD(K2,N)=FIELD(K2,N)-FIELD(K2,K1)
20    CONTINUE

C  SELECTION OF STARTING CONFIGURATIONS WITH NUCLEONS CLOSER THAN THE
C  HARD CORE RADIUS
      IF ((CORE.GT.0.).AND.(I.EQ.1)) CALL HCORE(*5,FIELD)
      Z=0.
      DO 25 K3=1,M
        Z=Z+FIELD(K3,N)**2
25    CONTINUE
      RSQU(I)=RSQU(I)+Z
      Y=(Y+Z)/B**2
      IF (Y.GT.200) Y=200
      R(I)=R(I)*EXP(-Y)

      RETURN
      END


      REAL FUNCTION RAN()

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE CALCULATES A RANDOM NUMBER.                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DOUBLE PRECISION X,XK,PM,A,C
      PARAMETER(A=14662125.D00,C=13136923.D00,PM=2.D00**48)

      COMMON/RG/X,XK

      XK=XK+.5D00
      X=A*X+DINT(C*XK)
      X=X-DINT(X/PM)*PM
      RAN=REAL(X/PM)
      RETURN
      END


      SUBROUTINE HCORE(*,FIELD)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS SUBROUTINE CHECKS WHETHER TWO PARTICLES ARE CLOSER THAN THE   C
C   GIVEN HARD CORE RADIUS. IF SO THE PROGRAM CONTINUES AT THE         C
C   ADDRESS SPECIFIED BY * IN THE CALLING PROGRAM, OTHERWISE STANDARD  C
C   RETURN TO THE CALLING PROGRAM                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(M=3,N=16)

      REAL FIELD(M,N)

      COMMON/HARDC/CORE

      J=0
      DO 20 K1=1,N-1
        DO 20 K2=K1+1,N
          SEPAR=0.
          DO 10 K3=1,M
            SEPAR=SEPAR+(FIELD(K3,K1)-FIELD(K3,K2))**2
10        CONTINUE
          IF (SEPAR.LT.CORE) RETURN 1
20    CONTINUE

800   RETURN
      END

      subroutine homo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS ROUTINE GENERATES CENTER OF MASS CORRELATED PARTICLES, WHICH  C
C   ARE HOMOGENEOUSLY DISTRIBUTED WITHIN A SPHERE, EXEPT FOR THIS      C
C   CORRELATION. THE RMS OF THE DISTRIBUTION IS THE SAME AS THE INPUT  C
C   RMS.                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      parameter (M=3,N=16)
      common/hardc/core
      real field(m,n)
      print*,'number of configurations ?'
      read*,numb
      print*,'physical rms             ?'
      read*,r0
      print*,' hard core radius        ?'
      read*,core
      core=core*core

C  CALCULATION OF THE RADIUS OF THE SPHERE; .64 IS THE SQUARE OF THE
C  EXTENSION OF THE NUCLEONS
      r0=sqrt(5.*(r0*r0-.64)/3.)
      pi=4.*atan(1.)
      write(3,1000)numb,1

C  LOOP
      do 100 k1=1,numb
1       field(1,n)=0.
        field(2,n)=0.
        field(3,n)=0.
        rsqu=0.
        do 50 k2=1,n-1
          phi=2.*pi*ran()
          ct=2.*(ran()-.5)
          st=sqrt(1.-ct*ct)
          r=r0*(ran()/3.)**(1./3.)
          field(1,k2)=r*st*cos(phi)
          field(2,k2)=r*st*sin(phi)
          field(3,k2)=r*ct
          field(1,n)=field(1,n)-field(1,k2)
          field(2,n)=field(2,n)-field(2,k2)
          field(3,n)=field(3,n)-field(3,k2)
          rsqu=rsqu+r*r
50      continue
        call hcore(*1,field)
        rsqu=rsqu+field(1,n)**2+field(2,n)**2+field(3,n)**2
        write(3,*)field
        write(3,1005)rsqu,1
100   continue
1000  format(2i8)
1005  format(e12.4,i8)
      return
      end


      SUBROUTINE OPTIC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE CALCULATES THE ELASTIC DIFFERENTIAL CROSS SECTION IN   C
C  THE OPTICAL POTENTIAL APPROXIMARION                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  IMPORTANT VARIABLES:                                                C
C     B     WIDTH PARAMETER OF THE GAUSSIAN DENSITY                    C
C     DR    DISTANCE BETWEEN MESHPOINTS IN INTEGRATION OF THE WAVE     C
C           FUNCTION                                                   C
C     F     ARRAY CONTAINING BESSEL AND NEUMANN FUNCTIONS              C
C     FTOT  SUM OF PARTIAL WAVE SCATTERING AMPLITUDES                  C
C     LMAX  NUMBER OF PARTIAL WAVES                                    C
C     L0    UPPER LIMIT FOR LMAX                                       C
C     NPTS  NUMBER OF STEPS IN THE INTEGRATION                         C
C     NSTEP ANGULAR STEP FOR DIFFERENTIAL CROSS SECTION                C
C     PL    LEGENDRE POLYNOMIAL                                        C
C     Q     MOMENTUM OF THE PROJECTILE                                 C
C     RHO0  NORMALIZATION OF THE POTENTIAL                             C
C     RK    MATCHING RADIUS TIMES Q                                    C
C     SLENG PROJECTILE NUCLEON SCATTERING LENGTH                       C
C     VEFF  (KINETIC ENERGY - POTENTIAL - ZENTRIFUGAL BARRIER)*DR**2   C
C     XD    LOGARITHMIC DERIVATIVE OF THE INTEGRATED WAVE FUNCTION     C
C     X/V   VALUES OF WAVE FUNCTION AND EFFECTIVE POTENTIAL AT AND     C
C            CLOSE TO THE MATCHING RADIUS                              C
C   AL1,B1,B2,PL0,PL1,Q1,V0,W,Y,X1,X2  AUXILIARY VARIABLES             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER(B=1.7D00,DR=4.D-2,N=16,NPTS=201,L0=20,NSTEP=9,kfm=10)

      COMPLEX X(5),V(5),Y(3),W(3),FL(0:L0),V0,RHO0s,rho0p,VEFF
      COMPLEX X1,X2,XD,IC,FTOT,ff,f0,f1,pleng(kfm),sleng(kfm)
      DOUBLE PRECISION F(0:L0+1,2)

      COMMON Q1,AL1,rho0s,rho0p,ddr,xi/BESS/F

1     CONTINUE
2     CONTINUE
      print*,' value of xi, xi=0 is without L.L.-correction '
      read*,xi
      xi=xi/3.d00

C  INPUT OF THE ANGULAR STEPS, MOMENTUM AND P-WAVE SCATTERING LENGTH
      READ(16,1045)KNUMB
      READ(16,*)SLENG
      READ(16,*)PLENG
      READ(16,*)Q
      READ(16,*)NTHET
      read(16,*)kf
      kf=min0(kf,kfmax)

      IC=CMPLX(0.,1.)
      PI=4.D00*DATAN(1.D00)
      Q1=Q*q
      RK=Q*(NPTS-1)*dr

      JUMP=4
4     WRITE(*,1100)L0
      READ(*,*,ERR=700)LMAX
      LMAX=IABS(LMAX)
      LMAX=MIN0(LMAX,L0)
      CALL BESSEL(LMAX,RK)


      do 500 j1=1,knumb
C
C   f0 defined with minus in potential, compensated here
C
         f0=-sleng(j1)
         f1=pleng(j1)
         ff=cmplx(12.d00*pi*n/q1)*f1
         rh=(3.d00*pi*b*b/4.d00)**(-1.5)
c  rho0s,rho0p sind die faktoren die die dichte multiplizieren sowie
c  die normierung der dichte
         RHO0s=cMPLX(4*PI*N*rh)*f0
         rho0p=cmplx(rh)*ff
         ddr=dr*dr

C  INTEGRATION OF THE RADIAL WAVE FUNCTION
      DO 200 K0=1,LMAX+1
        L=K0-1
        AL1=L*(L+1)*DR*DR
        Y(2)=0.
        Y(3)=CMPLX(DR**(L+1))
        V0=VEFF(2)
        W(2)=0.
        W(3)=Y(3)*(1.-V0/12.)
        DO 5 K1=3,NPTS+2
          W(1)=W(2)
          W(2)=W(3)
          Y(1)=Y(2)
          Y(2)=Y(3)
          W(3)=2.*W(2)-W(1)+V0*Y(2)
          V0=VEFF(K1)
          Y(3)=W(3)/(1.-V0/12.)
          IF (K1.GE.NPTS-2) THEN
            X(K1+3-NPTS)=Y(3)
            V(K1+3-NPTS)=V0
          END IF
5       CONTINUE

C  CALCULATION OF THE LOGARITHMIC DERIVATIVE OF THE INTEGRATED WAVE
C  FUNCTION AT THE MATCHING RADIUS R=DR*(NPTS-1)
        X1=(X(4)-X(2))*CMPLX(-8./21.)+(X(5)-X(1))*CMPLX(37./84.)
        X2=(V(4)*X(4)-V(2)*X(2))*CMPLX(-148./315.)+(V(5)*X(5)-V(1)*
     1       X(1))*CMPLX(-17./630.)
        XD=(X1+X2)/(CMPLX(DR)*X(3))

C  DETERMINATION OF THE PARTIAL WAVE SCATTERING AMPLITUDES FL(L)
       X1=CMPLX(F(L,1)*Q*(L+1))-CMPLX(RK*F(L+1,1)*Q)-XD*CMPLX(F(L,1)*RK)
       X2=CMPLX(-F(L,2)*Q*(L+1))+CMPLX(RK*Q*F(L+1,2))
       FL(L)=X1/(X2+XD*CMPLX(F(L,2)*RK)-IC*X1)
       FL(L)=CMPLX((2*L+1)/Q)*FL(L)

200   CONTINUE

C  STORING OF PARAMETERS AND PARTIAL WAVE SCATTERING AMPLITUDES
      WRITE(4,1020)f0,Q
      write(4,1021)f1,q
      write(15,1021)f1,q
      WRITE(4,1060)
      WRITE(4,1070)(K1-1,FL(K1-1),K1=1,LMAX+1)
      WRITE(4,1050)

C  CALCULATION OF THE DIFFERENTIAL CROSS SECTION FROM THE PARTIAL WAVE
C  SCATTERING AMPLITUDES USING A RECURSION RELATION FOR THE LEGENDRE
C  POLYNOMIALS

      IF (Q.GT.9.D-2) THEN
        DO 400 K1=1,181,NSTEP
          PL0=1.D00
          B1=DCOS(PI*(K1-1)/180)
          PL1=B1
          FTOT=FL(0)+FL(1)*CMPLX(PL1)
          DO 300 K2=2,LMAX
            B2=1.D00/K2
            PL=(2-B2)*B1*PL1+(B2-1)*PL0
            FTOT=FTOT+FL(K2)*CMPLX(PL)
            PL0=PL1
            PL1=PL
300       CONTINUE

C  ANGLE,SCATTERING AMPLITUDE AND CROSS SECTION ARE STORED
          WRITE(4,1000)K1-1,FTOT,10.*CABS(FTOT)**2
          WRITE(15,1005)j1,K1-1.,10.*CABS(FTOT)**2
400     CONTINUE
      ELSE
        FTOT=FL(0)
        WRITE(4,1000)0,FTOT,10.*CABS(FTOT)**2
      END IF

500   continue
      RETURN

700   PRINT*,'GIVE THE LAST INPUT ONCE MORE!  '
      GOTO (1,2),JUMP

1000  FORMAT(T9,I3,T20,E12.4,T34,E12.4,T53,E12.4)
1005  format(i3,2e12.4)
1045  format(i8)
1020  FORMAT(T5,'OPTICAL POTENTIAL CALCULATION:'/'    '
     1/T3,'SCATTERING LENGTH [FM]:',E12.4,' + I* ',E12.4/T3,'MOMENTUM [F
     2M**(-1)]:   ',E12.4/'  ')
1021  FORMAT(
     1T3,'SCATTERING LENGTH [FM]:',E12.4,' + I* ',E12.4/T3,'MOMENTUM [FM
     2**(-1)]:   ',E12.4/'  ')
1050  FORMAT('  '/'  '/T22,'SCATTERING AMPLITUDE [FM]'/T3,
     1'SCATTERING ANGLE',T22,'REAL PART',T33,'IMAGINARY PART',T52,'CROSS
     2 SECTION [mb]')
1060  FORMAT(T23,'(PARTIAL WAVE AMPLITUDES)*(2L+1)'/T19,'L',T34,'FL*(2L+
     11)')
1070  FORMAT(T17,I3,T23,E12.4,' + I* ',E12.4)
1100  FORMAT(T3,'MAXIMAL ANGULAR MOMENTUM (<=',I3,') ')

      END


      COMPLEX FUNCTION VEFF(K)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C THIS FUNCTION DETERMINES THE VALUE OF THE EFFECTIVE POTENTIAL        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (B=1.7D00,Dr=4.D-2)
      COMPLEX RHO0s,rho0p,cp,cdp,ar

      COMMON Q1,AL1,rho0s,rho0p,ddr,xi

      R=DR*(K-1)
      db=b*b
      R3=8*R/(db*3.D00)
      RHO1=0.D00
      IF (R3*r.LT.100) RHO1=DEXP(-R3*r/2.d00)
      ar=cmplx(rho1)*rho0p
      ar=ar/(1.-cmplx(xi)*ar)
      cp=ar*cmplx(r3)*(1.+ar*cmplx(xi))
      cdp=(cmplx(r3/r)-cmplx(r3*r3)-ar*cmplx(xi*(3.*r3*r3-r3/r))
     1    -cmplx(2.d00*(xi*r3)**2)*ar*ar)*ar
      VEFF=CMPLX(AL1/R**2.)-((cmplx(q1)-rho0s*cmplx(rho1)-0.5*cdp-
     1     cp/cmplx(r))/(1.-ar)+0.25*(cp/(1.-ar))**2)*cmplx(ddr)

      RETURN
      END


      SUBROUTINE BESSEL(LMAX,RK)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  IN THIS ROUTINE THE BESSEL AND NEUMANN FUNCTIONS ARE CALCULATED BY  C
C  BACKWARD AND FORWARD RECURSION RESPECTIVELY                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (L0=20)
      DOUBLE PRECISION F(0:L0+1,2),G(0:L0+30)

      COMMON/BESS/F

      C2=DSIN(RK)/RK
      F(0,2)=DCOS(RK)/RK
      F(1,2)=F(0,2)/RK+C2
      IF (LMAX.GE.1) THEN
        DO 10 K=2,LMAX+1
          F(K,2)=DBLE(2.*K-1.)*F(K-1,2)/RK-F(K-2,2)
10      CONTINUE
        G(LMAX+30)=0.D00
        G(LMAX+29)=1.D00
        DO 15 K2=1,LMAX+29
          G(LMAX+29-K2)=DBLE(2.*(LMAX-K2)+61)*G(LMAX+30-K2)/RK
          G(LMAX+29-K2)=G(LMAX+29-K2)-G(LMAX+31-K2)
15      CONTINUE
        DO 20 K3=1,LMAX+1
          F(K3,1)=C2*G(K3)/G(0)
20      CONTINUE
      ELSE
        F(1,1)=C2/RK-F(0,2)
      END IF
      F(0,1)=C2

      RETURN
      END


       subroutine f04adf(a,pnuml,co ,m1,m2,m3,b,m4,h,ifail)
C
C  Solve linear equations (call like NAGLIB f04adf)
C
C Maximum number of particles
       parameter(maxpnm=64)
C
C  Auxilliary variables
C
C  Counters etc.
       integer i,j,l,t,ai,aj,bi,bj,pnuml
       integer m1,m2,m3,m4,ifail
C  Greensfunction matrix
       complex a(pnuml,pnuml)
C  State vectors
       complex b(pnuml),ci(maxpnm),co(pnuml)
C Auxilliary variable for Gauss algorithm
       complex h(pnuml),hs

C
C  Third : Solve linear equations
C  Solution : b
C  Method : Bronstein, Taschenbuch der Mathematik
C           S.735 ( Gaussches Eliminationsverfahren )
C
c
c  save inhomogenous vector
c
          do 100 ai=1,pnuml
	     ci(ai)=co(ai) 
100       continue
c
c
c
C Transform A to triangular shape
          do 200 ai=1,(pnuml-1)
             do 230 aj=(ai+1),pnuml
                h(aj)=-a(aj,ai)/a(ai,ai)
                do 220 bi=1,pnuml
                  if (bi.eq.ai) then
                     a(aj,bi)=0.0
                  else
                     a(aj,bi)=a(aj,bi)+h(aj)*a(ai,bi)
                  endif
 220            continue
                ci(aj)=ci(aj)+h(aj)*ci(ai)
 230         continue
 200      continue
C Generate solution vector
          b(pnuml)=ci(pnuml)/a(pnuml,pnuml)
          do 300 ai=(pnuml-1),1,-1
             hs=0.0
             do 310 aj=1,(pnuml-ai)
                hs=hs+a(ai,ai+aj)*b(ai+aj)
 310         continue
          b(ai)=(ci(ai)-hs)/a(ai,ai)
 300      continue


	  end
