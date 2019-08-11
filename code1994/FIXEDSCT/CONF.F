      PROGRAM distri
c
c This program generates configuration of
c gaussian or homogenous type for
c a give width, number of particles and
c hardcore radius
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
C      M      DIMENSIONALITY OF THE SPACE                              C
C      N      NUMBER OF PARTICLES                                      C
C      B      WIDTH OF THE distribution function                       C
C   VARIABLES USED IN THE STOCHASTIC PART OF THE PROGRAM               C
C      NSTART  NUMBER OF INITIAL CONFIGURATIONS                        C
C      NKEEP   NUMBER OF CONFIGURATION KEPT AFTER EACH INITIAL         C
C              CONFIGURATION                                           C
C      X       INITIAL VALUE FOR RANDOM NUMBERS                        C
C      XD      COUNTER FOR RANDOM NUMBERS                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      parameter ( nmax = 50 , mmax = 3 )
      DOUBLE PRECISION X,XK
      real field(mmax,nmax)
 
      
      COMMON/RG/X,XK

      PI=4.*ATAN(1.)

C  INITIALIZATION OF THE RANDOM NUMBER GENERATOR
      X=PI*1.D11
      XK=0.D00
      DO 10 K1=1,20
        TT=RAN()
10    CONTINUE
      NAUS=1

      OPEN (3,FILE='CONFIG',STATUS='UNKNOWN')
      OPEN (16,FILE='cdata',STATUS='UNKNOWN')

      REWIND(3)
      REWIND(16)

      read(16,*) m
      read(16,*) n
      read(16,*) B
      read(16,*) core
      read(16,*) NSTART,NKEEP
      read(16,*) NT

      NSTART=NSTART*NKEEP
      core=core*core
      if (NT.eq.1) then
c
c generate gaussian distibution
c
        CALL GAUSS(field,m,n,b,core,NSTART)
      else
c
c generate homogenous distribution
c
        CALL homo(field,m,n,b,core,NSTART)
      endif
c
c calculate formfactor
c


      CLOSE (3,STATUS='KEEP')
      CLOSE (16,STATUS='KEEP')



1010  FORMAT(2I8)

      END



      SUBROUTINE GAUSS(field,m,n,b,core,numb)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS ROUTINE GENERATES CONFIGURATIONS WITH GAUSSIAN DISTRIBUTION   C
C   sinple product algoritm implemented                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      MULTI  MULTIPLICITY OF THE CONFIGURATION                        C
C      NCOUNT COUNTS FOR ACCEPTABLE CONFIGURATIONS                     C
C      NKEEP1 COUNTS THE CONFIGURATIONS STORED FOR AN INITIAL          C
C             CONFIGURATION                                            C
C      RAND   CONTAINS THE RANDOM NUMBERS CORRESPONDING TO THE         C
C             ACCEPTED AND TRIAL CONFIGURATION                         C
C      R      WEIGHT OF THE ACCEPTED AND TRIAL CONFIGURATION           C
C      RSQU   SUM OVER SQUARES OF PARTICLE DISTANCES                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real field(m,n)
      real b,core

      pi=4.*atan(1.)
      write(3,1000)numb,1

C  LOOP
      do 100 k1=1,numb
1       continue
        do 50 k2=1,n
          field(1,k2)=rgauss(b)
          field(2,k2)=rgauss(b)
          field(3,k2)=rgauss(b)
50      continue
        call hcore(*1,field,core,m,n)
        write(3,*)field
        write(3,1005)rsqu,1
100   continue
1000  format(2i8)
1005  format(e12.4,i8)
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


      SUBROUTINE HCORE(*,FIELD,core,m,n)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS SUBROUTINE CHECKS WHETHER TWO PARTICLES ARE CLOSER THAN THE   C
C   GIVEN HARD CORE RADIUS. IF SO THE PROGRAM CONTINUES AT THE         C
C   ADDRESS SPECIFIED BY * IN THE CALLING PROGRAM, OTHERWISE STANDARD  C
C   RETURN TO THE CALLING PROGRAM                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL FIELD(M,N)
      real core

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

      subroutine homo(field,m,n,b,core,numb)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   THIS ROUTINE GENERATES CENTER OF MASS CORRELATED PARTICLES, WHICH  C
C   ARE HOMOGENEOUSLY DISTRIBUTED WITHIN A SPHERE, EXEPT FOR THIS      C
C   CORRELATION. THE RMS OF THE DISTRIBUTION IS THE SAME AS THE INPUT  C
C   RMS.                                                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real field(m,n)
      real b,core

      pi=4.*atan(1.)
      write(3,1000)numb,1

C  LOOP
      do 100 k1=1,numb
1       continue
        do 50 k2=1,n
          phi=2.*pi*ran()
          ct=2.*(ran()-.5)
          st=sqrt(1.-ct*ct)
          r=1/b*ran()
          field(1,k2)=r*st*cos(phi)
          field(2,k2)=r*st*sin(phi)
          field(3,k2)=r*ct
          rsqu=rsqu+r*r
50      continue
        call hcore(*1,field,core,m,n)
        write(3,*)field
        write(3,1005)rsqu,1
100   continue
1000  format(2i8)
1005  format(e12.4,i8)
      return
      end

      real function rgauss(a)
C
C   Produces a gaussian
C   EXP(-x*x*a*a)
C   random number
C   (very ineffective!)
C
       rgauss=1/a*sqrt(-log(ran()))*cos(2*3.1415*ran())

       end
