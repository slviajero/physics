 	
        program main
          implicit none
          double precision xzero(100),xwei(100)
          integer n,i

           write(*,*) ' compute laguerre weights & points: '
           write(*,*) ' give N:'
           read (*,*)   n
           call lag_point(n,xzero,xwei)
           write (*,*) 'i,  xzero,   xwei*exp(xzero)'
           do 1, i=1,n
                write(*,*) i, xzero(i), xwei(i)

1          continue
           end


	SUBROUTINE  LAG_POINT(N,XZERO,XWEI)
************************************************************
* written by      :  THOMAS KRAUS
*
* purpose         :  computes N zeroes and weights of the
*                    Nth Laguerre-polynomial - the zeroes
*                    are stored in XZERO, while the weights
*                    are put into XWEI.
*                    The weights are multiplied by exp(XZERO).
*
* last change     :  3/24/92
*
*************************************************************
      IMPLICIT NONE
**************Lokale Variablen: *****************************
	SAVE
      INTEGER    N
      DOUBLE PRECISION  XZERO(100),XWEI(100)
      DOUBLE PRECISION  EPSILON,XACC,XLOW,XUP,FAC
      PARAMETER(EPSILON=.0001)
      INTEGER    J,ICOUNT,M
      DOUBLE PRECISION  ALAGUERRE,RTSAFE
      EXTERNAL   RTSAFE,ALAGUERRE,FUN_LAG

C     check whether N>1:
      IF(N.LE.1) THEN
           WRITE(*,*) ' N TOO SMALL IN LAG_POINT !!'
           STOP
      ENDIF
      XZERO(1)=.5857864376
	XZERO(2)=3.4142135623
      ICOUNT=2
C     Now begin with recurrence relations:
C     XACC is the accuracy for the root
      XACC=0.00030
      DO 1 ICOUNT=2,N-1
C          M defines the Laguerre polynomial for which roots are just
C          searched
           M=ICOUNT+1
           IF(ICOUNT.EQ.N-1) XACC=1.D-10
C               Add new root:
                XLOW=XZERO(ICOUNT)+EPSILON
                XUP =2.*XLOW
                XZERO(ICOUNT+1)=RTSAFE(FUN_LAG,M,XLOW,XUP,XACC)
                DO 2,J=ICOUNT,2,-1
                     XLOW= XZERO(J-1)+EPSILON
                     XUP = XZERO(J)  -EPSILON
                     XZERO(J)=RTSAFE(FUN_LAG,M,XLOW,XUP,XACC)
 2              CONTINUE
		    XLOW = 0.0000001D0  
		    XUP  = XZERO(1)-EPSILON
                XZERO(1)=RTSAFE(FUN_LAG,M,XLOW,XUP,XACC)
 1    CONTINUE
C     Finally compute weights (cp Abramowitz p.890):
C     First calculate N! :
      FAC=1.D0
      DO 6,J=1,N
           FAC=DFLOAT(FAC)*DFLOAT(J)
 6    CONTINUE
      M=N+1
      DO 5,J=1,N
           XWEI(J)=DEXP(XZERO(J))*XZERO(J)/DFLOAT(M)**2
     +              /(ALAGUERRE(M,XZERO(J)))**2
 5    CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION ALAGUERRE(L,X)
************************************************************
*
* program         :  DOUBLE PRECISION FUNCTION ALAGUERRRE
*
* written by      :  Thomas Kraus
*
* purpose         :  creation of Laguerre-Polynomials by
*                    recurrence relations
*
* last change     :  3/23/92
*
*************************************************************
      IMPLICIT NONE
	SAVE
****** arguments: *******************************************
      DOUBLE PRECISION X
      INTEGER L
**************local variables: *****************************
      DOUBLE PRECISION P0,P1,PDUMMY
      INTEGER I

**** see whether range of L and X is correct: ****************
      IF ((L.LT.0).OR.(X.LT.0.0)) THEN
           WRITE(*,*) ' X OR L OUT OF RANGE IN ALAGUERRE!!'
           STOP
      ENDIF
      IF (L.EQ.0) THEN
          ALAGUERRE = 1.D0
      ELSE
           P0=1.D0
           P1=1.D0-X
           DO 10, I=2,L
                PDUMMY=P1
                P1= 2*P1 - P0 - ((X+1.D0)*P1 -P0)/I
                P0= PDUMMY
10         CONTINUE
           ALAGUERRE = P1
      ENDIF
      RETURN
      END
	DOUBLE PRECISION FUNCTION RTSAFE(FUNCD,L,X1,X2,XACC)
C     from NUMERICAL RECIPES p258
C     Using a combination of Newton_Raphson and bisection, find
c     the root of a function bracketed between X1 and X2. The 
C     root, returned as the function value RTSAFE , will be re-
C     fined until its accuracy is known within +-XACC. FUNCD is
C     a user supplied subroutine which returns both the function
C     value and the first derivative of the function.
C     Attention !!!!!
C     in this implementation the subroutine FUNCD must be of the 
C     form FUNCD(L,X,VALUE,DERIVATIVE) with L being an integer
C     specifying the function. (the other arguments are DOUBLE
C	PRECISION).
C***************************************************************
	IMPLICIT NONE
	SAVE
C     arguments:
	DOUBLE PRECISION X1,X2,XACC
	INTEGER L,J
C     maximum number of allowed iterations:
	INTEGER MAXIT
	PARAMETER(MAXIT=500)
	DOUBLE PRECISION FL ,FH,DF,XH,DX,DXOLD,F,XL,TEMP
	EXTERNAL FUNCD

	CALL FUNCD(L,X1,FL,DF)
	CALL FUNCD(L,X2,FH,DF)
	IF(FL*FH.GE.0.) PAUSE ' ROOT MUST BE BRACKETED IN RTSAFE !'
	IF(FL.LT.0.) THEN
		XL=X1
		XH=X2
	ELSE
		XH=X1
		XL=X2
	ENDIF
	RTSAFE=.5*(X1+X2)
	DXOLD=ABS(X2-X1)
	DX=DXOLD
	CALL FUNCD(L,RTSAFE,F,DF)
	DO 11 J=1,MAXIT
		IF(((RTSAFE-XH)*DF-F)*((RTSAFE-XL)*DF-F).GE.0.
     +           .OR.ABS(2.*F).GT.ABS(DXOLD*DF)) THEN
			DXOLD=DX
			DX=0.5*(XH-XL)
			RTSAFE=XL+DX
			IF(XL.EQ.RTSAFE) RETURN
      	ELSE 
			DXOLD=DX
			DX=F/DF
			TEMP=RTSAFE 
			RTSAFE=RTSAFE-DX
			IF(TEMP.EQ.RTSAFE) RETURN
		ENDIF
		IF (ABS(DX).LT.XACC) RETURN
		CALL FUNCD(L,RTSAFE,F,DF)
		IF (F.LT.0.) THEN
			XL=RTSAFE
		ELSE
			XH=RTSAFE
	      ENDIF
 11	CONTINUE 
	PAUSE 'RTSAFE EXCEEDING MAXIMUM ITERATIONS !'
	END
	SUBROUTINE FUN_LAG(L,X,Y,DV)
C     computes the Lth Laguerre-Polynomial at argument X
C     value is returned in Y, the first derivative is stored
C     in DV
C************************************************************
	IMPLICIT NONE
	SAVE
	DOUBLE PRECISION X,Y,DV
	INTEGER L,M
	DOUBLE PRECISION ALAGUERRE
	EXTERNAL ALAGUERRE

	
C     Check range of X and L:
	IF((X.LE.0).OR.(L.LT.0)) THEN
		WRITE(*,*) ' VALUES OUT OF RANGE IN FUN_LAG !!'
		STOP
      ENDIF

	Y=ALAGUERRE(L,X)
	IF(X.LT.1.D-11) THEN
		DV=0.D0
	ELSE 
		M=L-1
		DV=L/X*(ALAGUERRE(L,X)-ALAGUERRE(M,X))
	ENDIF
	RETURN
	END

