      PROGRAM srdbl

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
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PARAMETER USED IN THE WHOLE PROGRAM:                               C
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
c
c 30.07.1992 : - Konfigurationserzeugung
c                fuer 2,3,4 Streuzentren eingebaut
c              - n Abhaengigkeit der stochastischen
c                Rechnung eingebaut
c              - s,p und sp Programme zusammenkopiert
c              - Optisches Potential entfernt
c              - Hardcore auf n Teilchen umgestellt
c 01.08.1992   - Probleme mit conv1 und conv2 unter
c                Profortran
c                geloest durch schreiben einer eigenen
c                sinus funktion
c 02.08.1992   - Problem der numerischen Instabilitaet
c                in der p-Welle
c                (1) Ursache gefunden Im pleng << Im h_2
c                (2) Programm geaendert pleng wird nicht
c                    mehr rausdividiert
c                (3) Umschreiben auf Double
c
c
c 12.10.1992   - Programmversion mit Konfigurations
c                erzeugung waehrend der Rechnung erstellt
c                dazu sp und p Welle entfernt
c                bei Konfigurationserzeugung alles bis
c                auf homogene Dichte entfernt
c 26.10.1992   - variabler Radius in homogene Dichte 
c                programmiert
c
      parameter(M=3,Nm=200)

      double precision x,xk
      real pi
      real r0
      real core
      integer nstart,nkeep,n
      real tt
      integer k1,naus
c
c global variables for matrices
c
      real*8 field(M,Nm)
      complex*16 a(Nm,Nm)
c
c common blocks for random numbers (rg)
c and general parameters (para)
c
      COMMON/PARA/PI,r0,NSTART,NKEEP
      COMMON/RG/X,XK
      common/hardc/core
c
c open confdata file and read configuration data information
c
      open(21,file='confdata')
c number of scattering centers
      read(21,*) n
      if ((n.gt.200).or.(n.lt.2)) then
         print *,'Illegal number of scattering centers'
         stop
      endif
c number of configurations
      read(21,*) nstart
      nkeep=1
c hardcore radius
      read(21,*) core
      core=core*core
c r0 read from file
      read(21,*) r0
      close(21)
c
c pi loaded
c
      PI=4.*ATAN(1.)
c
c parameters for homogenous distribution
c
C  CALCULATION OF THE RADIUS OF THE SPHERE; .64 IS THE SQUARE OF THE
C  EXTENSION OF THE NUCLEONS
c  this was done for oxygen. As we are  not interested in a constant
c  density r0 is not scaled but the correct value for oxygen is used

      r0=sqrt(5.*(r0*r0-.64)/3.)
c
C  INITIALIZATION OF THE RANDOM NUMBER GENERATOR
c
      X=PI*1.D11
      XK=0.D00
      DO 10 K1=1,20
         TT=RAN()
10    CONTINUE
      NAUS=1

      OPEN (1,FILE='form',STATUS='UNKNOWN')
      OPEN (2,FILE='stoch',STATUS='UNKNOWN')
      OPEN (3,FILE='config',STATUS='UNKNOWN')
      OPEN (15,FILE='flucti',STATUS='UNKNOWN')
      OPEN (16,FILE='f1data',STATUS='UNKNOWN')
      OPEN (17,FILE='forwar',STATUS='UNKNOWN')
      open(18,file='cross',status='unknown')
      open(20,file='report')

      REWIND(1)
      REWIND(2)
      REWIND(3)
      REWIND(15)
      REWIND(16)
      REWIND(17)

      CALL EDCS_s(n,field,a)

300   CLOSE (1,STATUS='KEEP')
      CLOSE (2,STATUS='KEEP')
      CLOSE (3,STATUS='KEEP')
      CLOSE (4,STATUS='KEEP')
      CLOSE (15,STATUS='KEEP')
      CLOSE (16,STATUS='KEEP')
      CLOSE (17,STATUS='KEEP')

800   STOP


1010  FORMAT(2I8)

      END



      SUBROUTINE EDCS_s(n,field,am)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE CALCULATES THE SCATTERING AMPLITUDE AT 0 DEGREES IF    C
C  THE MOMENTUM OF THE INCOMING PARTICLE IS SMALLER THAN 0.2 INVERSE   C
C  FERMI; OTHERWISE DIFFERENTIAL CROSS SECTION. IN ADDITION, THE       C
C  NUMBER OF DIFFERENT CONFIGURATIONS, THE AVERAGE ACCEPTANCE          C
C  PROPABILITY, AND THE RMS-RADIUS IS DETERMINED.                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PARAMETER:                                                         C
C      KAMAX   NUMBER OF VALUES OF THE SCATTERING ANGLE WHERE THE      C
C              DIFFERENTIAL CROSS SECTION IS CALCULATED                C
C      KMAX    UPPER LIMIT FOR THE NUMBER OF VALUES OF THE INCIDENT    C
C              MOMENTUM                                                C
C      Q0      UPPER LIMIT FOR CALCULATION OF SCATTERING AMPLITUDE IN  C
C              FORWARD DIRECTION ONLY                                  C
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

      PARAMETER(M=3,nm=200)
      parameter(KAMAX=61,KMAX=10,Q0=.0001,kfm=10)

c
c variables which homo loads with
c configuration data
c
      REAL*8 FIELD(M,N)
      real*8 rsqu
c
c internal data of s-wave routine
c
      real*8 RIJ(Nm,Nm),QANGLE(2,KAMAX)
      COMPLEX*16 IC,AM(N,N),PHI(Nm),FAKT,FKKP(KMAX,KAMAX),
     1        PHI1(Nm,KAMAX),sleng0,pleng(kmax),sleng(kmax)
c
c parameters from main program
c
      real pi
      real r0
      integer nstart,nkeep
      COMMON/PARA/PI,r0,NSTART,NKEEP

C  INPUT OF THE ANGULAR STEPS, MOMENTUM AND P-WAVE SCATTERING LENGTH
      READ(16,1045)KNUMB
      read(16,*)sleng
      READ(16,*)PLENG
      NTHET=180
      READ(16,*)Q
      READ(16,*)NTHET
      NTHET=MAX0(4,NTHET)
      mthet=nthet
      kangle=181
      CONTINUE
      read(16,*)kf
      kf=min0(kf,kfm)


      IC=CMPLX(0.,1.)
      SUM=0.
      RMS=0.

C  CALCULATION OF UNIT VECTORS SPECIFYING THE DIRECTION OF THE
C  OUTGOING PROJECTILE
      K10=0
      DO 8 K1=1,181,NTHET
	K10=K10+1
	QANGLE(1,K10)=SIN(PI*(K1-1)/180.)
	QANGLE(2,K10)=COS(PI*(K1-1)/180.)
	DO 8 K2=1,KNUMB
	  FKKP(K2,K10)=CMPLX(0.)
8     CONTINUE
c
c central loop:
c    nstart is the number of configurations
c    nkeep and multi are one
c    this construction is kept only for historical
c    reasons  ( see MUSH )
c
       multi=1
       DO 200 K1=1,NSTART*NKEEP
            call homo2(field,rsqu,n)
c
C  BEGINNING OF THE CALCULATION OF THE SCATTERING AMPLITUDE FOR ONE
C  CONFIGURATION
c
	SUM=SUM+MULTI
	RMS=RMS+RSQU*MULTI
c
C  DETERMINATION OF THE PARTICLE DISTANCES
c
	DO 20 K2=1,N-1
	  DO 20 K3=K2+1,N
	    R=0.
	    DO 10 K4=1,M
	      R=R+(FIELD(K4,K2)-FIELD(K4,K3))**2
10          CONTINUE
	    RIJ(K2,K3)=R**.5
20      CONTINUE
c
C   START OF THE CALCULATION OF SCATTERING AMPLITUDES
c
	DO 200 K2=1,KNUMB
c
c scattering length taken from f1data
c
	  sleng0=sleng(k2)
c
C   CALCULATION OF THE INHOMOGENEOUS VECTOR AND THE OUTGOING WAVES
c
	  DO 30 K3=1,N
	    PHI(K3)=EXP(IC*Q*FIELD(M,K3))
	    K10=0
	    DO 30 K4=1,KANGLE,MTHET
	      K10=K10+1
	      PHI1(K3,K10)=EXP(-IC*Q*(FIELD(M,K3)*QANGLE(2,K10)+
     1                       FIELD(M-1,K3)*QANGLE(1,K10)))
30      CONTINUE
c
C   CALCULATION OF THE MATRIX ELEMENTS
c
	  DO 40 K3=1,N-1
	    AM(K3,K3)=1.
	    DO 40 K4=K3+1,N
	      AM(K3,K4)=-SLENG0*EXP(IC*Q*RIJ(K3,K4))/RIJ(K3,K4)
	      AM(K4,K3)=AM(K3,K4)
40        CONTINUE
	  AM(N,N)=CMPLX(1.)
c
C   SOLUTION OF THE SYSTEM OF MULTIPLE SCATTERING EQUATIONS
c   attention this routine works only for symmetric matrices
c
	  DO 70 K3=1,N-1
	    DO 70 K4=K3+1,N
	      FAKT=AM(K3,K4)/AM(K3,K3)
	      AM(K3,K4)=FAKT
	      DO 70 K5=K3+1,N
		AM(K5,K4)=AM(K5,K4)-AM(K5,K3)*FAKT
70        CONTINUE
	  DO 90 K4=2,N
	    DO 90 K5=1,K4-1
	      PHI(K4)=PHI(K4)-AM(K5,K4)*PHI(K5)
90        CONTINUE
	  PHI(N)=PHI(N)/AM(N,N)
	  DO 110 K4=1,N-1
	    K3=N-K4
	    DO 100 K5=K3+1,N
	      PHI(K3)=PHI(K3)-AM(K5,K3)*PHI(K5)
100         CONTINUE
	    PHI(K3)=PHI(K3)/AM(K3,K3)
110       CONTINUE
c
C   AVERAGING OF THE INDIVIDUAL SCATTERING AMPLITUDES
c
	  DO 200 K6=1,N
	    K10=0
	    DO 200 K5=1,KANGLE,MTHET
	      K10=K10+1
	      FKKP(K2,K10)=FKKP(K2,K10)+MULTI*PHI(K6)*PHI1(K6,K10)
200   CONTINUE
c
C  STORING OF SCATTERING ANGLES, AVERAGED SCATTERING AMPLITUDES AND
C  DIFFERENTIAL CROSS SECTION
c
      DO 300 K1=1,KNUMB
	sleng0=sleng(k1)
	WRITE(2,1020)n,SLENG0
	WRITE(2,1030)Q
	IF (K1.EQ.1) WRITE(2,1050)
	K10=0
	write(17,1100) 98,abs(sleng0),real(fkkp(k1,1)*sleng0)/sum
	write(17,1100) 99,abs(sleng0),imag(fkkp(k1,1)*sleng0)/sum
	write(17,1100) 96,q,real(fkkp(k1,1)*sleng0)/sum
	write(17,1100) 97,q,imag(fkkp(k1,1)*sleng0)/sum
	gint=0.
	DO 301 K2=1,KANGLE,MTHET
	  K10=K10+1
	  WRITE(2,1040)K2-1,SLENG0*FKKP(K1,K10)/SUM,10.*(ABS(FKKP(K1,
     1          K10)*SLENG0)/SUM)**2
	  gint=gint+10.*abs(fkkp(k1,k10)*sleng0/sum)**2
     1           *sin((k2-1)/180.*pi)*mthet/180.*pi
301     continue
	write(2,*) 'Elastic cross section   [mb]',gint*2*pi
	write(2,*) 'Total cross section     [mb]',
     1               40*pi/q*imag(fkkp(k1,1)*sleng0/sum)
	write(2,*) 'Inelastic cross section [mb]',
     1               40*pi/q*imag(fkkp(k1,1)*sleng0/sum)-gint*2*pi
	write(2,*) 'Elementary cross section [mb]',
     1               40*pi/q*imag(sleng0)

	write(18,1100)1,q,gint*2*pi
	write(18,1100)2,q,40*pi/q*imag(fkkp(k1,1)*sleng0/sum)
	write(18,1100)3,q,40*pi/q*imag(fkkp(k1,1)*sleng0/sum)-gint*2*pi
	write(18,1100)4,q,40*pi/q*imag(sleng0)
300   CONTINUE
c
C   STORING OF INFORMATION CONCERNING THE CONFIGURATIONS AND THE RMS; THE
C   .64 BEING ADDED IS THE SQUARE OF THE NUCLEON RADIUS
c
      WRITE(2,1060)NSTART*NKEEP,NSTART*NKEEP/SUM,(RMS/(SUM*N)+.64)**.5
      write(2,*) 'radius of homogenous density: ',r0*(1./3.)**(1./3.)
800   RETURN


1005  FORMAT(E12.4)
1010  FORMAT(E12.4,I8)
1020  FORMAT(T5,'STOCHASTIC CALCULATION ',i3,' CENTERS :'/'    '
     1/T3,'SCATTERING AMPLITUDE [FM]:',E12.4,' + I* ',E12.4,'  ')
1030  FORMAT(T3,'MOMENTUM [FM**(-1)]:   ',E12.4/'  ')
1045  FORMAT(I3)
1050  FORMAT(T22,'SCATTERING AMPLITUDE [FM]'/T3,
     1'SCATTERING ANGLE',T22,'REAL PART',T33,'IMAGINARY PART',T52,'CROSS
     2 SECTION [mb]')
1060  FORMAT(' '/T5,'NUMBER OF DIFFERENT CONFIGURATIONS',T60,I8/T5,
     1'AVERAGE ACCEPTANCE PROPABILITY',T60,E12.4/T5,'RMS-RADIUS',T60,
     2E12.4/'  ')
1040  FORMAT(T8,I4,T20,2E12.4,T52,E12.4)
1090  FORMAT(T5,' # OF VALUES FOR THE INCIDENT MOMENTUM K ( <= ',I2,')
     1 ')
1100  FORMAT(I3,2E12.4)
1110  FORMAT(I3,4E12.4)

      END




      SUBROUTINE HCORE(*,FIELD,n)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS SUBROUTINE CHECKS WHETHER TWO PARTICLES ARE CLOSER THAN THE   C
C   GIVEN HARD CORE RADIUS. IF SO THE PROGRAM CONTINUES AT THE         C
C   ADDRESS SPECIFIED BY * IN THE CALLING PROGRAM, OTHERWISE STANDARD  C
C   RETURN TO THE CALLING PROGRAM                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      PARAMETER(M=3)

      REAL*8 FIELD(M,N)

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

 


      subroutine homo2(field,rsqu,n)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS ROUTINE GENERATES CENTER OF MASS CORRELATED PARTICLES if      C
C   N < 4 else noncorellated particles, WHICH                          C
C   ARE HOMOGENEOUSLY DISTRIBUTED WITHIN A SPHERE, EXEPT FOR THIS      C
C   CORRELATION. THE RMS OF THE DISTRIBUTION IS approximatly           C
C   THE SAME AS THE INPUT                                              C
C   RMS.                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c  - input : r0, core, n
c  - output: field,rsqu
c  every time this subroutine is called a new
c  homogenous configuration is generated and stored
c  into field
c

      parameter (M=3)
      real pi
      integer nstart,nkeep
      real ran
      real*8 field(m,n)
      real*8 rsqu
      common/hardc/core
      common/para/pi,r0,nstart,nkeep
c
2         continue
	  rsqu=0.
	  do 150 k2=1,n
	    phi=2.*pi*ran()
	    ct=2.*(ran()-.5)
	    st=sqrt(1.-ct*ct)
	    r=r0*(ran()/3.)**(1./3.)
	    field(1,k2)=r*st*cos(phi)
	    field(2,k2)=r*st*sin(phi)
	    field(3,k2)=r*ct
	    rsqu=rsqu+r*r
150       continue
	  call hcore(*2,field,n)
c
c now field contains the configuration
c and rsqu then square radius
c
1000    format(2i8)
1005    format(e12.4,i8)
	return
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

