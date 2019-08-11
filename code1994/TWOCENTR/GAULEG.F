	program main
       implicit none
        integer i,n
        double precision xlow,xup,x(100),w(100)

        write(*,*) 'weights & points for gauss-legendre:'
        write(*,*) 'give N:'
        read(*,*) n
	  xlow=-1.d0
        xup=1.d0
        call gauleg (xlow,xup,x,w,n)
         write(*,*) 'i,  xzero,   xweight'
        do 1,i=1,n
           write(*,*) i,x(i),w(i)
 1       continue
         end

********************************************************
*  									 *
*       3 JUNE 1992                                    *
*                                                      *
*       SUBROUTINE GAULEG(X1,X2,X,W,N)                 *
*                                                      *
*       Given the lower and upper limits of integration*
*       X1 and X2, and given N, this routine returns   *
*       arrays X and W of length N, containing the     *
*	  abscissas and weights of the Gauss_Legendre    *
*       N-point quadrature formula.                    *
*                                                      *
*       ( Numerical Recipes p. 125 )                   *
*                                                      *
********************************************************

	SUBROUTINE GAULEG(X1,X2,X,W,N)

	IMPLICIT REAL*8 (A-Z)
	SAVE
	INTEGER I,J,M,N
	PARAMETER (EPS=3.D-14)
	DIMENSION X(N),W(N)


	M=(N+1)/2
	XM=0.5D0*(X1+X2)
	XL=0.5D0*(X2-X1)
	DO 12 I=1,M
		Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
 
 1		CONTINUE
			P1=1.D0
			P2=0.D0
			DO 11 J=1,N
				P3=P2
				P2=P1
				P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
 11			CONTINUE

			PP=N*(Z*P1-P2)/(Z*Z-1.D0)
			Z1=Z
			Z=Z1-P1/PP

		IF(ABS(Z-Z1).GT.EPS) GO TO 1
		X(I)=XM-XL*Z
		X(N+1-I)=XM+XL*Z
		W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
		W(N+1-I)=W(I)
 12	CONTINUE
	RETURN
	END
