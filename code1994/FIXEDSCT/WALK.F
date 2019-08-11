       program walk
       integer nconf
       parameter(nconf=20000)
       real x(2*nconf,3),y(2*nconf,3)
       integer nx,ny
       integer i,j,k,k1
       real v0,vpot,r,b
       real m,dt
       integer nt
       real tt,pi
       real fakt
       real ct,st,phi
       real ran,rgau
       external ran,rgau
       
       double precision xr,xk
       COMMON/RG/Xr,XK
       
       open(10,file='out')
       pi =4.*atan(1.)
c
C  INITIALIZATION OF THE RANDOM NUMBER GENERATOR
      Xr=PI*1.D11
      XK=0.D00
      DO 10 K1=1,20
        TT=RAN()
10    CONTINUE
c
c Eingabe der Parameter 
      print *,'Zahl der Schritte :'
      read *,nt
      print *,'dt :'
      read *,dt
      print *,'Weite des Potentials :'
      read *,b
      print *,'Tiefe des Potentials :'
      read *,v0
      print *,'Masse des Projectils :'
      read *,m
c
c Anfangskonfiguration
c wichtige Bemerkung:      
c das Potential das hier gesampelt wit enthaelt einen 
c Normierungsfaktor
c von der Normierung der gaussschen Zufallszahlen.
c das muss spaeter korrigiert werden.
c
c gauss:
       do 20 j=1,nconf
          do 20 k=1,3
             x(j,k)=rgau()*b
20    continue
c kasten:
c       do 20 j=1,nconf
c          phi=2.*pi*ran()
c          ct=2.*(ran()-.5)
c          st=sqrt(1.-ct*ct)
c          r=b*(ran())**(1./3.)
c          x(j,1)=r*st*cos(phi)
c          x(j,2)=r*st*sin(phi)
c          x(j,3)=r*ct
c20     continue     
      nx=nconf
      fakt=1.
c
c Propagation
      do 100 i=1,nt
c
c Propagation: sample Gausssche Zufallszahlen
c mit Weite m/(dt)
       do 110 j=1,nx
         do 110 k=1,3
           y(j,k)=sqrt(dt/m)*rgau()+x(j,k)
110    continue
       ny=nx
c dieses Porgram geht nur fuer repulsive Potentiale
c Annihilation:
       nx=1
       do 120 j=1,ny
             r=y(j,1)**2+y(j,2)**2+y(j,3)**2
c gauss:
             vpot = v0*exp(-r/(b**2)/2.)
c kasten:
c            if (sqrt(r).le.b) vpot=v0
             tt   = exp(-dt*vpot)
c
c Modifikation fuer attraktive Potentiale
c 
99           if (tt.gt.1.) then
                x(nx,1)=y(j,1)
                x(nx,2)=y(j,2)
                x(nx,3)=y(j,3)
                nx=nx+1 
		tt=tt-1
		goto 99
             endif 
             if (ran().lt.tt) then
                x(nx,1)=y(j,1)
                x(nx,2)=y(j,2)
                x(nx,3)=y(j,3)
                nx=nx+1 
             endif
120    continue
       nx=nx-1
c
c Kontrolle ob nicht zu viele vernichtet
c oder erzeugt
c wenn ja, verdoppelung
       if (nx.lt.nconf/2) then
          do 130 j=1,nx
             do 130 k=1,3
               x(j+nx,k)=x(j,k)
130       continue
          fakt=fakt*2
          nx=nx*2
       endif 
c oder reduzierung
c im Mittel werden die Haelfte gestrichen
c
       if (nx.gt.(3*nconf)/2) then
          ny=0
	  do 131 j=1,nx 
             if (ran().gt.(.5)) then
		ny=ny+1
		do 132 k=1,3
		   x(ny,k)=x(j,k)
132             continue
	     endif 
131       continue
	  fakt=(fakt*ny)/nx
	  nx=ny
        endif


c       print *,dt*i,nx/fakt
c
c berechne Ergebniss
       if (mod(i,10).eq.0) then
         tt=0.
          do 140 j=1,nx
            r=x(j,1)**2+x(j,2)**2+x(j,3)**2
c  gauss:
              tt = tt+exp(-r/(b**2)/2.)
c kasten:
c             if (sqrt(r).le.b) tt = tt + 1.
140      continue
         tt=tt/(nconf*fakt)/sqrt(2*pi)*sqrt(m)*(i*dt)**(3./2.)
c gauss:
          tt=v0**2*(2.*pi)**(3./2.)*b**3*tt
c kasten:
c        tt=v0**2*(4*pi/3.*b**3)*tt
         print *,i*dt,(abs(tt))*4*pi 
         write(10,*) i*dt,(abs(tt))*4*pi
       endif
100    continue
       close(10)              
       end

      real function rgau()
      real pi,x1,x2
      real ran
      external ran
      pi=4.*atan(1.)
      x1=ran()
      x2=ran()
      rgau=sqrt(-2.*alog(x1))*cos(2*pi*x2)
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


