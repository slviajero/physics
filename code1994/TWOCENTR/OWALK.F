      PROGRAM o4
	call optic
      END

      SUBROUTINE OPTIC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(DR=4.D-2,N=4,NPTS=201,L0=20,NSTEP=9,kfm=10)
      COMPLEX X(5),V(5),Y(3),W(3),FL(0:L0),V0,VEFF
      COMPLEX X1,X2,XD,IC,FTOT
      DOUBLE PRECISION F(0:L0+1,2)
      real b,vo0,m,al1
      COMMON b,vo0,m,al1/BESS/F

      print *,'Weite des Potentials :'
      read *,b
      print *,'Tiefe des Potentials :'
      read *,vo0
      print *,'Masse des Projektils :'
      read *,m

      q=0.01

      IC=CMPLX(0.,1.)
      PI=4.D00*DATAN(1.D00)
      Q1=Q*q
      RK=Q*(NPTS-1)*dr

      JUMP=4
      LMAX=0
      LMAX=IABS(LMAX)
      LMAX=MIN0(LMAX,L0)
      CALL BESSEL(LMAX,RK)

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
      WRITE(*,1060)
      WRITE(*,1070)(K1-1,FL(K1-1),K1=1,LMAX+1)
      WRITE(*,1050)

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
	  WRITE(*,1000)K1-1,FTOT,10.*CABS(FTOT)**2
400     CONTINUE
      ELSE
	FTOT=FL(0)
	WRITE(*,1000)0,FTOT,10.*CABS(FTOT)**2
      END IF
      
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (Dr=4.D-2)
      real b,v0,m
      real al1,r
      COMMON b,v0,m,al1
      r=dr*(k-1)
c gauss
      VEFF=CMPLX(AL1/R**2.+2.*m*v0*exp(-(r/b)**2/2)*dr*dr)
c kasten      
c       if (r.gt.b) then   
c           veff=cmplx(0.)
c       else
c           veff=cmplx(2*m*v0*dr*dr)
c       endif
c       veff=veff+cmplx(al1/r**2)
      
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




