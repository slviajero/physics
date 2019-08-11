
      PROGRAM o16   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  ORGANISATION OF THE PROGRAM:                                        C
C    SUBROUTINE OPTIC CALCULATES ELASTIC DIFFERENTIAL CROSS SECTION IN C
C                     OPTICAL MODEL APPROXIMATION                      C
C    COMPLEX FUNCTION VEFF  PROVIDES THE VALUES OF THE EFFECTIVE       C
C                     POTENTIAL FOR INTEGRATING THE SCHROEDINGER EQ.   C
C    SUBROUTINE BESSEL DETERMINES THE VALUES OF BESSEL AND NEUMANN     C
C                      FUNCTIONS                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PARAMETER USED IN THE WHOLE PROGRAM:                               C                          C
C      N      NUMBER OF PARTICLES                                      C
C      B      WIDTH OF THE HARMONIC OSCILLATOR WAVE FUNCTION           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   O16:
C   Authors: D.Stoll and S. Lenz
C   Version: 26.11.91
C
C  This optical potential program is based on Dieter Stolls
C  MUSH. It can calculate the scattering amplitude, the partial
C  waves and the diff. cross section at momentum k from
C  a given optical potential VEFFA.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C version message
C
      character*(80) version
      version='O4 optical potential program, version 1.0, 26.11.91'
C
C files for input and output
C
      OPEN (4,FILE='optical',STATUS='UNKNOWN')
      OPEN (15,FILE='flucti',STATUS='UNKNOWN')
      OPEN (16,FILE='f1data',STATUS='UNKNOWN')
      OPEN (17,FILE='forwar',STATUS='UNKNOWN')

      REWIND(4)
      REWIND(15)
      REWIND(16)
      REWIND(17)


        print *,version
100     continue
           CALL OPTIC
           print *,'Enter 1 for new calculation'
           read *,i
           if (i.ne.1) goto 999
        goto 100
999     continue


      CLOSE (4,STATUS='KEEP')
      CLOSE (15,STATUS='KEEP')
      CLOSE (16,STATUS='KEEP')
      CLOSE (17,STATUS='KEEP')



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
C     SLENG PROJECTILE NUCLEON s-wave SCATTERING AMPLITUDE             C
C     PLENG PROJECTILE NUCLEON p-wave SCATTERING AMPLITUDE             C
C     VEFF  (KINETIC ENERGY - POTENTIAL - ZENTRIFUGAL BARRIER)*DR**2   C
C     XD    LOGARITHMIC DERIVATIVE OF THE INTEGRATED WAVE FUNCTION     C
C     X/V   VALUES OF WAVE FUNCTION AND EFFECTIVE POTENTIAL AT AND     C
C            CLOSE TO THE MATCHING RADIUS                              C
C   AL1,B1,B2,PL0,PL1,Q1,V0,W,Y,X1,X2  AUXILIARY VARIABLES             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER(B=1.7D00,DR=4.D-2,N=16,NPTS=201,L0=20,kfm=10)

      integer nstep
      real sreact,selast
      COMPLEX X(5),V(5),Y(3),W(3),FL(0:L0),V0,RHO0s,rho0p,VEFFA
      COMPLEX X1,X2,XD,IC,FTOT,fll,fv,f0,f1,pleng(kfm),sleng(kfm)
      DOUBLE PRECISION F(0:L0+1,2)

      COMMON Q1,AL1,rho0s,rho0p,ddr,xi/BESS/F


C       print *,' Plonummer: '
C       read *,np

      JUMP=1
1     continue
2     continue
C
C  Lorentz-lorentz correction
C  not implemented now
C
      print*,' value of xi, xi=0 is without L.L.-correction '
      read*,xi
      xi=xi/3.d00
C
C  INPUT OF THE ANGULAR STEPS, MOMENTUM AND P-WAVE SCATTERING LENGTH
C
      READ(16,*)KNUMB
      READ(16,*)SLENG
      READ(16,*)PLENG
      READ(16,*)Q
      READ(16,*)NSTEP
      read(16,*)kf
      kf=min0(kf,kfm)

      IC=CMPLX(0.,1.)
      PI=4.D00*DATAN(1.D00)
      Q1=Q*Q
      RK=Q*(NPTS-1)*dr

      JUMP=4
4     WRITE(*,1100)L0
      READ(*,*,ERR=700)LMAX
      LMAX=IABS(LMAX)
      LMAX=MIN0(LMAX,L0)
      CALL BESSEL(LMAX,RK)

      do 500 j1=1,knumb
      f0=SLENG(j1)
      f1=pleng(j1)
C
C   fll is the s-wave scattering length  (including all factors
C         that come with it)
C   fv is the p-wave scattering volume
C   rh is the normalisation factor of the gaussian
C
      fll=cMPLX(4*PI*(N-1.d00))*f0
      fv=-cmplx(12.d00*pi*(n-1.d00)/q1)*f1
      rh=0.25*(1.0D00/(b*b*pi))**(1.5)
C
c  rho0s,rho0p  are the factors of the density in veff
c  they are calculated once and handed down in a common block
C
      RHO0s=cMPLX(rh)*fll
      rho0p=cmplx(rh)*fv
      ddr=dr*dr

C  INTEGRATION OF THE RADIAL WAVE FUNCTION
      DO 200 K0=1,LMAX+1
        L=K0-1
C
C angular momentum term al1
C
        AL1=L*(L+1)*ddr
        Y(2)=0.
        Y(3)=CMPLX(DR**(L+1))
        V0=VEFFA(2)
        W(2)=0.
C
C  ????
C
        W(3)=Y(3)*(1.-V0/12.)
        DO 5 K1=3,NPTS+2
          W(1)=W(2)
          W(2)=W(3)
          Y(1)=Y(2)
          Y(2)=Y(3)
          W(3)=2.*W(2)-W(1)+V0*Y(2)
          V0=VEFFA(K1)
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
      WRITE(4,1020)f0,q
      write(4,1021)f1,q
      write(15,1021)f1,q
      WRITE(4,1060)
      WRITE(4,1070)(K1-1,FL(K1-1),K1=1,LMAX+1)
      WRITE(4,1050)

C  CALCULATION OF THE DIFFERENTIAL CROSS SECTION FROM THE PARTIAL WAVE
C  SCATTERING AMPLITUDES USING A RECURSION RELATION FOR THE LEGENDRE
C  POLYNOMIALS

      IF (Q.GT.9.D-2) THEN
c
c   scattering amplitude and differential cross section
C   are clculated here
c   selast contains elastic crosssecton
c   sreact contains imaginary part of forw. scatt. ampl times 
c   4 pi over k
c

        sreact=0.
	selast=0.
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
          if (k1.eq.1) then 
	     sreact=4.*pi/q*aimag(ftot)
          endif
          selast=10.*CABS(FTOT)**2*sin((k1-1)/180.*pi)*nstep/180.*pi+
     a           selast
	  WRITE(4,1000)K1-1,FTOT,10.*CABS(FTOT)**2
          write(1,1200) q*(2.*(1.001-cos((k1-1)/180.*pi)))**.5,ftot
          WRITE(15,1005)j1,K1-1.,10.*CABS(FTOT)**2
          if ((k1-1).eq.0) then
           write(17,2000) np,cabs(f1),aimag(ftot)
           write(17,2000) np+1,cabs(f1),real(ftot)
2000       format(i3,2e12.4)
          endif
400     CONTINUE
      ELSE
        FTOT=FL(0)*cmplx(n/(n-1.d00))
        sreact=4.*pi/q*aimag(ftot)
	WRITE(4,1000)0,FTOT,10.*CABS(FTOT)**2
      END IF
      write(4,401) 10.*sreact
      write(4,402) selast*2*pi
401   format('reaction crosssection : ',e12.4)
402   format('elastic  crosssection : ',e12.4)



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
1200  FORMAT(d12.4,2e12.4)
      END



      COMPLEX FUNCTION VEFFA(K)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C THIS FUNCTION DETERMINES THE VALUE OF THE EFFECTIVE POTENTIAL        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  This is an alternative effective potential subroutine
C  it allows for easy change of the used density
C  The potential is supossed to be in Kisslinger shape
C  (->Helgert).
C
C  the user has to give:
C    the s-wave rhos density as a function of r,q1 and f0
C    the p-wave rhop density as a function of r,q1 and f1  (and xi)
C        and its first two derivatives rhop1 and rhop2
C   ( the densities must be already multiplied by
C     the scattering amplitudes, particle numbers
C     see Klaus Helgerts thesis )
C----------------------------------------------------------------------
C
C  b is the width of the gaussian
C  dr ist the step length
C
      implicit logical (a-z)
      DOUBLE PRECISION B,DR,db
      PARAMETER (B=1.7D00,Dr=4.D-2)
      COMPLEX RHO0s,rho0p,cp,cdp,ar,frk
      COMPLEX rhos,rhop,rhop1,rhop2,ff1,ff2,part1,part2
      COMPLEX gamma,h1,h2
C-----------------------------------------------------------------------
C  parameters are handed down by this common block
C  q1 is the momentum
C  rho0s and rho0p are the scattering amplitudes
C    times the normalisation of the density etc
C  ddr is dr*2
C  al1 is the angular momentum factor, it contains
C  already the ddr factor !!!
C
C
      DOUBLE PRECISION q1,al1,ddr,xi,r,r3,rho1
      COMMON Q1,AL1,rho0s,rho0p,ddr,xi
      real h4
C-----------------------------------------------------------------------
C
C r : distance from origin
      R=DR*(K-1)
C parmeter for the gaussian
      db=b*b
      R3=R/db
C the density rho1 is set to 0 outside some range
      RHO1=0.D00
      IF (R3*r.LT.100) RHO1=DEXP(-R3*r)
C
C here the densities are defined
C in a very ineffective way
      h4=1. 
C s-wave density
      rhos=rho0s*cmplx(rho1*(1.D00+h4*2.*r3*r))
C p-wave (without LL correction)
      rhop=rho0p*cmplx(rho1*(1.D00+h4*2.*r3*r))
      rhop1=-2.*cmplx(r3)*rhop+h4*4.*rho0p*cmplx(rho1*r3)
      rhop2=-cmplx(2./db)*rhop-cmplx(r3*2.)*rhop1+
     1       h4*(cmplx(4./db*rho1)*rho0p-cmplx(8.*r3*r3*rho1)*rho0p)
C
C p-wave LL corrected
C the LL correction replaces the p-wave density rhop
C (attention ! including all factors that come with it)
C by rhop/(1+xi*rhop) where c=-1/3
C
c     gamma=1.
c      /(1.+cmplx(xi)*rhop)
c      h1=rhop*gamma
c      h2=rhop1*gamma**2
c      rhop2=rhop2*gamma**2-2.*gamma*cmplx(xi)*rhop1*h2
c      rhop1=h2
c      rhop=h1
C
C effective potential
C  some functions:
      frk=1.+rhop
C
C part 1 : determination of the s-wave part
C
      part1=+(cmplx(q1)+rhos)/frk
C
C part2 : determination of the kissliger part
C
      part2=-.25*(cmplx(2./r)+rhop1/frk)**2-
     1       .5*(cmplx(-2./(r*r))+rhop2/frk-(rhop1/frk)**2)
C
C the effective potential needs an additional minus
C
      veffa=-(part1+part2)*cmplx(ddr)+cmplx(al1/(r*r))

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
