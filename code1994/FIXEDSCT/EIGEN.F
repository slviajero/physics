        program eigen

	real*8 a(100,100)
	integer n
	real*8 e(100)
	real*8 r(100)
	integer i,j,ifail

        open(10,file='res')

	do 50 n=4,96,4

	do 110 i=1,n
	  do 110 j=1,n
             a(i,j)=0.
110     continue

	do 120 i=1,n
	     a(i,i+1)=1.
	     a(i+1,i)=1.
120     continue

        ifail=0

	call f02aaf(a,100,n,r,e,ifail)

	print *,' emax :', r(n)-2.

	write(10,*) n,abs(r(n)-2.)

50      continue

	end
