       program sk16
C
C Calculates scattering amplitudes for small k
C
       parameter(n=16,fmax=10,b=1.41)
       integer fnum,nthet,nav
       integer i,j,k
       complex sleng(fmax),pleng(fmax)
       complex f1,f1ren
       real q,ct

       open(10,file='STOCHASIC')
       open(11,file='forwar')
       open(12,file='F1DATA')


       rewind(10)
       rewind(11)
       rewind(12)

       read(12,*) fnum
       read(12,*) sleng
       read(12,*) pleng
       read(12,*) q
       read(12,*) nthet
       read(12,*) nav


       do 1000 i=1,fnum
         f1=pleng(i)
         f1ren=f1/(1-(0.0,1.0)*4.*sqr(2.)*f1*(n-1)/b)
         write(11,2000) 1,cabs(f1),3*n*f1ren
         do 1000 j=0,180,nthet
            ct=cos(pi/180*j)
            write(12,2000) i-1   ,j,cabs(3*n*f1ren*ct)**2*10
            write(12,2000) i-1+10,j,real(3*n*f1ren*ct)
            write(12,2000) i-1+20,j,aimag(3*n*f1ren*ct)
1000   continue

2000   format(i3,2e12.4)

       end