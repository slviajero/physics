       subroutine f04adf(a,pnuml,co ,m1,m2,m3,b,m4,h,ifail)
C
C  Solve linear equations (call like NAGLIB f04adf) 
C
C Maximum number of particles
       parameter(maxpnm=100)
C
C  Auxilliary variables
C
C  Counters etc.
       integer i,j,l,t,ai,aj,bi,bj
       integer m1,m2,m3,m4,ifail
C  Greensfunction matrix
       complex a(maxpnm,maxpnm)
C  State vectors
       complex b(maxpnm),ci(maxpnm),co(maxpnm)
C Auxilliary variable for Gauss algorithm
       complex h(maxpnm),hs

       integer pnuml

C
C  Third : Solve linear equations
C  Solution : b
C  Method : Bronstein, Taschenbuch der Mathematik
C           S.735 ( Gaussches Eliminationsverfahren )
C
c
c  save inhomogenous vector
c
          do 100 ai=1,pnuml
	     ci(ai)=co(ai) 
100       continue
c
c
c
C Transform A to triangular shape
          do 200 ai=1,(pnuml-1)
             do 230 aj=(ai+1),pnuml
                h(aj)=-a(aj,ai)/a(ai,ai)
                do 220 bi=1,pnuml
                  if (bi.eq.ai) then
                     a(aj,bi)=0.0
                  else
                     a(aj,bi)=a(aj,bi)+h(aj)*a(ai,bi)
                  endif
 220            continue
                ci(aj)=ci(aj)+h(aj)*ci(ai)
 230         continue
 200      continue
C Generate solution vector
          b(pnuml)=ci(pnuml)/a(pnuml,pnuml)
          do 300 ai=(pnuml-1),1,-1
             hs=0.0
             do 310 aj=1,(pnuml-ai)
                hs=hs+a(ai,ai+aj)*b(ai+aj)
 310         continue
          b(ai)=(ci(ai)-hs)/a(ai,ai)
 300      continue


	  end
