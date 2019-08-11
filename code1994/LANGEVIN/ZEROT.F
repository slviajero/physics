         program test
          
	 real x,dx,dt,dw
         real ds

	 print *,' dt = '
	 read *,dt
         ds = sqrt(2.*dt)

	 x=rgau()

	 do i=1,1000000
	    dx=f(x)*dt+ds*rgau()
	    x=x+dx
	    if (mod(i,1000).eq.0) print *,i,x
	    if (x.lt.1) print *,'*',i,x
         end do 

         end
         
	 real function f(x)

         real d,x,dt 

	 d=.1*(x**2-2.)*2.*x/((1.-x)*(1.+x))
         f=d
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


      real function rgau()
      real x1,x2
      real ran
      external ran
      real pi
      common/const/pi
      x1=ran()
      x2=ran()
      rgau=sqrt(-2.*alog(x1))*cos(2.*pi*x2)
      end
