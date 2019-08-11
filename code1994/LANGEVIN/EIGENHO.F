        program eigen

	real*8 a(500,500)
	integer n
	real*8 e(500)
	real*8 r(500)
	integer i,j,ifail
        real lambda

        open(10,file='res')

        b=100
	om2=.01

	do 50 n=10,400,10

        eps=.5
	alpha=eps/b**2

	do 110 i=1,n
	  do 110 j=1,n
             a(i,j)=0.
110     continue

	do 120 i=1,n
	     a(i,i+1)=1.
             a(i,i)  =-om2*eps**2
	     a(i+1,i)=1.
120     continue


        a(1,1)=1.-alpha-om2*eps**2/2.
        a(n,n)=1.-alpha-om2*eps**2/2.


        ifail=0

	call f02aaf(a,500,n,r,e,ifail)

	lambda=(r(n)-2.)/(eps)
	print *,' emax :',n*eps,lambda,abs(1/lambda)
	write(10,*) n*eps,lambda,abs(1/lambda)

50      continue

	end
