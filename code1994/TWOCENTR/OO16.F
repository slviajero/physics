       program optx
       implicit none
c parameter:
       integer npts
       parameter(npts=800)
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
       complex vopt
       external vopt
       complex psip,dpsip
       complex k1,k2,k3
       real rp
       real dx
       complex qp,xe,pp
C
C triviales:
       integer i
       real pi
       real x
       complex ic
C Standarddefinitionen:
C
C eingabedaten:
C m:    Projetilmasse
       real m,b
       complex v0
       complex f0
       common /param/m,b,v0,pi

       pi = 4.*atan(1.)
       ic = (0. , 1.)
C
C liest die Parameter vom Eingabefile
C
       b=1.7
       print *,'Streulaenge'
       read *,f0
       m=1.
       v0=4*pi*15*f0     
c
c Bestimmung der Reichweite des Potentials
c und der Schrittweite
c
	rp=b*4.
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

       print *,'Streulaenge: ',pp*(16./15.)

       end


       complex function vopt(x)
       real x
C
C lokales:
       real pi
       real m,b
       complex v0
       
       common /param/m,b,v0,pi
       
       if ((x/b)**2.lt.40.) then       
	   xe=exp(-(x/b)**2)
       else
	   xe=0.
       endif
       vopt=.25*(1./(b**2*pi))**(3./2.)*v0*xe*(1+2.*(x/b)**2)
       end




