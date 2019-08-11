       program optx
       implicit none
c
c Autor: 
C Datum:
C Version: 0.0
C-----------------------------------------------------------
C Berechnet die Steulaenge (bei k=0) fuewr beliebige Potentiale
C Diese Version benutzt die Runge Kutta Methode aus Abramowitz
C S897
C-----------------------------------------------------------
C
C Das Potential 
c
C v: elementares Potential
c
       complex*16 v
       external v
c
c effektives Potential
c
       complex*16 veff
       external veff
c
c optisches Potential
c
       complex*16 vop
       external vop
C
C eingabedaten:
C

       real*8 m,b,v0,omega,af0
       integer a
       real*8 pi
       common /param/m,b,v0,pi,omega,af0,a

       real*8 alpha

       complex*16 pp,ppe

       complex*16 runge
       external runge

       pi = 4.*atan(1.)

c
c Ausgabefiles
c 
       open(10,file="wave")
       open(11,file="oow-result")

10     continue

C
C liest die Parameter vom Eingabefile
C
       print *,'Weite des Potentials:'
       read *,b
       if (b.eq.0) goto 9999
       print *,'Staerke des Potentials:'
       read *,v0
       print *,'Masse '
       read *,m

c
c Korrektur mit harmonischen Oszillator Wellenfunktionen
c Annahme: beide Massen sind gleich!
c

       a=0
       omega=0

       print *,'Omega '
       read *,omega

       print *,'A '
       read *,a

       alpha=2*b**2*m*omega

       print *,"Alpha=", alpha

c
c Bestimmung der Reichweite des Potentials
c und der Schrittweite
c
c
c Berechnung der logarithmischen Ableitungen:
c
c
c Berechnung der Streulaenge:
c
       pp = runge(v,b) 
       ppe=pp

       print *,'Streulaenge elementar: ',pp
       print *,'Querschnitt elementar: ',4*pi*abs(pp)**2

       pp = runge(veff,b*sqrt((alpha+1)/alpha)) 

       print *,'Streulaenge effekiv: ',pp
       print *,'Querschnitt effektiv: ',4*pi*abs(pp)**2

c
c Vorbereitung für optisches Potential
c

       b=b/2
       m=m*2

       pp=runge(v,b)
       af0=dreal(pp)

       print *,"a fuer vopt", af0

       m=m/2
       b=b*2

       pp = runge(vop,sqrt(1/(2*m*omega)))

       print *,'Streulaenge optical: ',pp
       print *,'Querschnitt optical: ',4*pi*abs(pp)**2

       goto 10

9999   continue

       close(10)
       close(11)

       end

c
c Das elementare Potential
c

       complex*16 function v(x)
       real*8 x
C
C lokales:
C

       real*8 m,b,v0,omega,af0
       integer a
       real*8 pi
       common /param/m,b,v0,pi,omega,af0,a

       complex *16 xe

       if ((x/b)**2.lt.80.) then       
	   xe=exp(-(x/b)**2/2)
       else
	   xe=0.
       endif
       v=2.*v0*m*xe
       end

c
c Das effektive Potential
c

       complex*16 function veff(x)
       real*8 x
C
C lokales:
C
       real*8 pi
       real*8 m,b,v0,omega
       real*8 af0
       complex *16 xe
       integer a
c
c die Veraenderung des elementaren Potentials
c
       real*8 alpha
       real*8 bs,v0s

       
       common /param/m,b,v0,pi,omega,af0,a

       alpha=2*b**2*m*omega


       bs=b*sqrt((alpha+1)/alpha)
       v0s=v0*a*sqrt(alpha/(alpha+1))**3

       if ((x/bs)**2.lt.80.) then       
	   xe=exp(-(x/bs)**2/2)
       else
	   xe=0.
       endif
       veff=2.*v0s*m*xe
       end
c
c Das optische Potential
c

       complex*16 function vop(x)
       real*8 x
C
C lokales:
C
       real*8 pi
       real*8 m,b,v0,omega
       real*8 af0
       complex *16 xe
       integer a
c
c Parameter des opt. Potential
c
       real*8 bs,v0s

       
       common /param/m,b,v0,pi,omega,af0,a


       bs=1/sqrt(2*omega*m)
       v0s=-af0*(a-1)*2*sqrt(pi*omega/m)

       if ((x/bs)**2.lt.80.) then       
	   xe=exp(-(x/bs)**2/2)
       else
	   xe=0.
       endif
       vop=2.*v0s*m*xe
       end


C
C Die Runge Kutta Integration
C

      complex*16 function runge(vopt,b)
      implicit none
      complex*16 vopt
      external vopt
c
c typische Reichweite des Potentials
c
       real*8 b

c parameter:
       integer npts
       parameter(npts=4000)
c
C npts: Zahl der Integrationspunkte
c
C psip  : Wellenfunktion 
C dpsi  : Ableitung dazu
C k1,k2,k3 : Runge Kutta Hilfsgroessen
C rp  : Reichweite des Potentials       
C dx  : Schrittweite
C qp  : logarithmische Ableitung am Rand       
c
       complex*16 psip,dpsip
       complex*16 k1,k2,k3
       real*8 rp
       real*8 dx
       complex*16 qp

c
c triviales
c
      real*8 x
      integer i
      logical debug
      parameter (debug=.false.)
c
c integrationsdaten
c

	rp=b*6.
	dx=rp/npts
c       
c Intergrationsroutine:
c Anfangsbedingungen:
c
       psip  = (0.,0.)
       dpsip = sqrt(-vopt(0.0D00))
c
c Integration:
       do 100 i=1,npts
	  x=i*dx

          if (debug.and.(mod(i,100).eq.0)) then
		print *,i,x,psip,dpsip,vopt(x)
          endif
		
c
c Runge Kutta Hilfsgroessen bestimmen
	  k1=dx*vopt(x)*psip
	  k2=dx*vopt(x+dx/2)*(psip+dx/2*dpsip+dx/8*k1)
	  k3=dx*vopt(x+dx)*(psip+dx*dpsip+dx/2*k2)
c
c Interation
c
	  psip = psip+dx*(dpsip+(k1+2.*k2)/6.)
	  dpsip = dpsip+k1/6.+2./3.*k2+k3/6.
	  
100    continue

       qp=dpsip/psip
c
c Berechnung der Streulaenge:
c
       runge = 1/qp - rp 

      end
