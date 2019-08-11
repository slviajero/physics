       program optx
       implicit none
c parameter:
       integer npts
       parameter(npts=4000)
C
C npts: Zahl der Integrationspunkte
C-----------------------------------------------------------
C Autor:   Stefan Lenz
C Datum:
C Version: 0.0
C-----------------------------------------------------------
C Berechnet die Steulaenge (bei k=0) fuewr beliebige Potentiale
C Diese Version benutzt die Runge Kutta Methode aus Abramowitz
C S897
C zusaetzlich wurde eine A/(A-1) Korrektur eingefuegt
C-----------------------------------------------------------
C
C lokales:
C vopt: optisches Potential
C psip  : Wellenfunktion 
C dpsi  : Ableitung dazu
C k1,k2,k3 : Runge Kutta Hilfsgroessen
C rp  : Reichweite des Potentials       
C dx  : Schrittweite
C qp  : logarithmische Ableitung am Rand       
       complex*16 vopt
       external vopt
       complex*16 psip,dpsip
       complex*16 k1,k2,k3
       real*8 rp
       real*8 dx
       complex*16 qp,pp
C
C triviales:
       integer i
       real*8 pi
       real*8 x
       complex*16 ic
C Standarddefinitionen:
C
C eingabedaten:
C m:    Projetilmasse
       real*8 m,b,v0

       common /param/m,b,v0,pi

       pi = 4.*atan(1.)
       ic = (0. , 1.)
C
C liest die Parameter vom Eingabefile
C
       print *,'Weite des Potentials:'
       read *,b
       print *,'Staerke des Potentials:'
       read *,v0
       print *,'Masse '
       read *,m
      
c
c Bestimmung der Reichweite des Potentials
c und der Schrittweite
c
	rp=b*6.
	dx=rp/npts
	print *,'dx = ',dx
c       
c Intergrationsroutine:
c Anfangsbedingungen:
       psip  = (0.,0.)
       dpsip = sqrt(-vopt(0.))
c Integration:
       do 100 i=1,npts
	  x=i*dx
c
c Runge Kutta Hilfsgroessen bestimmen
	  k1=dx*vopt(x)*psip
	  k2=dx*vopt(x+dx/2)*(psip+dx/2*dpsip+dx/8*k1)
	  k3=dx*vopt(x+dx)*(psip+dx*dpsip+dx/2*k2)
c
c Interation
	  psip = psip+dx*(dpsip+(k1+2.*k2)/6.)
	 dpsip = dpsip+k1/6.+2./3.*k2+k3/6.
100    continue
c
c Berechnung der logarithmischen Ableitungen:
c
       qp=dpsip/psip
c
c Berechnung der Streulaenge:
c
       pp = 1/qp - rp 

       print *,'Streulaenge: ',pp
       print *,'Querschnitt: ',4*pi*real(pp)**2

       end


       complex*16 function vopt(x)
       real*8 x
C
C lokales:
       real*8 pi
       real*8 m,b,v0
       complex *16 xe
       
       common /param/m,b,v0,pi
       
       if ((x/b)**2.lt.80.) then       
	   xe=exp(-(x/b)**2/2)
       else
	   xe=0.
       endif
       vopt=2.*v0*m*xe
       end




