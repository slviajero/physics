       program gborn
	 

c 
c this program calculates the interacting
c averaged greensfunction g(t,x,y) at x=y=0 in
c first order born approximation
c -------------------------------------------------
c ATTENTION : - 1/ (4 pi ) is not included in
c              the energy dependent greensfunction
c              1 / ( 2 pi ) is not included into 
c              time fourier transform
c--------------------------------------------------
c for N scattering centers which are distributed
c according to a gaussian distribution
c
c      exp( - x*x * b*b )
c
c the scattering amplitude is unitarized with
c
c    a(E) = (exp(2 I delta ) - 1 ) /(2 I * sqr(E) )
c
c Input data (read from file 10 = ginp )
c     N
c     b
c     delta
c     tres,btres,tstep
c     
c Output data (written to file 11 = gout )
c (in PLOT Format)
c plot 1  : Re G
c plot 2  : Im G
c plot 3  : abs(G)**2
c plot 11 : Re G0
c plot 12 : Im G0
c plot 13 : abs(G0)**2
c
c
c input data:
       integer N
       real b,tres,btres,tstep
       complex delta
c variables:
       real t
       complex g,g0
c
c parameters
       complex IC
       real pi
       IC = (0.0 , 1.0)
       pi = 3.1415
c
c read input data:
       open(10,FILE='ginp')
       open(11,FILE='gout')
       read (10,*) N
       read (10,*) b
       read (10,*) delta
       read (10,*) tres,btres,tstep            
c
c calculate greensfunction
       do 10 i=1,50
	   t=tres*(btres)**i
           g0 = (2*pi) * (IC*sqrt(pi)) * 1/csqrt(IC*t)**3
           g  = g0+
     A            2*N*b*b*(exp(2*IC*delta)-1)/sqrt(IC*t+1/(b*b))*
     B            ( sqrt(pi)/(2*IC)
     C              +log((sqrt(1+IC*t*b*b)+1)/(sqrt(1+IC*t*b*b)-1))
     D            )
            write(11,1000) 1,t,real(g)
            write(11,1000) 2,t,Imag(g)
            write(11,1000) 3,t,abs(g)**2
 
            write(11,1000) 11,t,real(g0)
            write(11,1000) 12,t,Imag(g0)
            write(11,1000) 13,t,abs(g0)**2
10       continue 

1000    format(i3,2e12.4)

        end
 



       
