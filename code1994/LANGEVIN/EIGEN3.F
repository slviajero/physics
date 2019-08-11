        program eigen

	real*8 a(400,400)
	integer n
	real*8 e(400)
	real*8 r(400)
	integer i,j,ifail

        open(10,file='res')

        n=200

	do b=.1,10.,.1

	do 110 i=1,n
	  do 110 j=1,n
             a(i,j)=0.
110     continue

	do 120 i=1,n
	     a(i,i+1)=1.
	     a(i+1,i)=1.
120     continue
   
	a(1,1)=1.-.5/b**2
	a(n,n)=1.-.5/b**2

        ifail=0

	call f02aaf(a,400,n,r,e,ifail)

	print *,' emax :',b, r(n)-2.

	write(10,*) b,abs(r(n)-2.)

        end do 

	end
