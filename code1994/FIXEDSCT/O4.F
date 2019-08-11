      PROGRAM P4new

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  ORGANISATION OF THE PROGRAM:                                        C
C    SUBROUTINE EDCS  CALCULATES ELASTIC DIFFERENTIAL CROSS SECTION    C
C    SUBROUTINE FORM  CALCULATES THE FORM FACTOR OF THE DISTRIBUTION   C
C    SUBROUTINE GAUSS GENERATES A GAUSSIAN DENSITY DISTRIBUTION        C
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
      OPEN (17,FILE='forwar',STATUS='UNKNOWN')
      open (20,file='life' )

      REWIND(1)
      REWIND(2)
      REWIND(3)
      REWIND(4)
      REWIND(15)
      REWIND(16)
      REWIND(17)
      rewind(20)

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
      ELSEIF (NANT.EQ.3) THEN
        CALL OPTIC
      ELSEIF (NANT.EQ.2) THEN
        CALL homo
        naus=2
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
        CALL EDCS_p
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
1001  FORMAT(T8,' 3 : OPTICAL POTENTIAL CALCULATION'
     1/T6,' ELSE  END THE PROGRAM'/'   ')
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


      SUBROUTINE EDCS_p

      END



      SUBROUTINE FORM

      END



      SUBROUTINE GAUSS
      end



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

      PARAMETER(B=1.41D00,DR=4.D-2,N=4,NPTS=201,L0=20,NSTEP=9,kfm=10)

      COMPLEX X(5),V(5),Y(3),W(3),FL(0:L0),V0,RHO0s,rho0p,VEFF
      COMPLEX X1,X2,XD,IC,FTOT,ff,f0,f1,pleng(kfm),sleng(kfm)
      DOUBLE PRECISION F(0:L0+1,2)

      COMMON Q1,AL1,rho0s,rho0p,ddr,xi/BESS/F


       np=100
      xi=-1.0/3.d00

C  INPUT OF THE ANGULAR STEPS, MOMENTUM AND P-WAVE SCATTERING LENGTH
      READ(16,1045)KNUMB
      read(16,*)sleng
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
      f1=pleng(j1)
      f0=sleng(j1) 
      ff=cmplx(12.d00*pi*(n-1.d00)/q1)*f1
      rh=(3.d00*pi*b*b/4.d00)**(-1.5)
c  rho0s,rho0p sind die faktoren die die dichte multiplizieren sowie
c  die normierung der dichte
      RHO0s=cMPLX(4*PI*(N-1.d00)*rh)*f0
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

      IF (Q.GT.1.001D-1) THEN
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
          ftot=ftot*cmplx(n/(n-1.d00))
          WRITE(4,1000)K1-1,FTOT,10.*CABS(FTOT)**2
          WRITE(15,1005)j1,K1-1.,10.*CABS(FTOT)**2
          if ((k1-1).eq.0) then
           write(17,2000) np,cabs(f1),aimag(ftot)
           write(17,2000) np+1,cabs(f1),real(ftot)
2000       format(i3,2e12.4)
          endif
400     CONTINUE
      ELSE
        FTOT=FL(0)*cmplx(n/(n-1.d00))
        WRITE(4,1000)0,FTOT,10.*CABS(FTOT)**2
        write(20,*) sngl(q),atan(aimag(ftot)/real(ftot))
      END IF

500   continue
      RETURN
700   continue
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
      PARAMETER (B=1.41D00,Dr=4.D-2)
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
       parameter(maxpnm=100)
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


