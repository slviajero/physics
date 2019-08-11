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
C Loest das optische Potential fuer einen harmonischen
C Oszillator im Impulsraum.
C Diese Version benutzt die Runge Kutta Methode aus Abramowitz
C S897
C zusaetzlich wurde eine A/(A-1) Korrektur eingefuegt
C-----------------------------------------------------------
C
C lokales:
C vopt: optisches Potential
C psim  : Wellenfunktion Paritaet -
C dpsim : Ableitung dazu
C psip  : Wellenfunktion Paritaet +
C dpsi  : Ableitung dazu
C k1,k2,k3 : Runge Kutta Hilfsgroessen
C v0  : Potential - Hilfsvariable       
C rp  : Reichweite des Potentials       
C dx  : Schrittweite
C qp  : logarithmische Ableitung der P+ Welle am Rand       
C qm  : logarithmische Ableitung der P- Welle am Rand
C t   : Tramsmissionsamplitude       
C r   : Reflektionsamplitude
       complex vopt
       external vopt
       complex psim,dpsim
       complex psip,dpsip
       complex k1,k2,k3
       complex v0
       real rp
       real dx
       complex qp,qm
       complex r,t
C
C triviales:
       integer i
       real pi
       real x
       complex ic
       complex pp,pm
       character*(100) line
C Standarddefinitionen:
C
C eingabedaten:
C mt:    Targetmasse
C mp:    Projetilmasse
C hb:    Planck-const
C omega: Oszillatorstaerke
C pin:   einlaufender Impuls
C tau:   elementare WW-Staerke
C nt:    Zahl der Targetzustaende
C np:    Zahl der Impulszustaende
       real mp,mt
       real hb
       real omega
       real pin
       complex tau
       integer nt,np

       common /param/mp,mt,hb,omega,pin,tau,nt,np,pi

       pi = 4.*atan(1.)
       ic = (0. , 1.)
C
C liest die Parameter vom Eingabefile
C
       open(10,file='parin')
       read(10,*) mt,mp
       read(10,*) hb
       read(10,*) omega
       read(10,*) pin
       read(10,*) tau
       read(10,*) nt,np
       close(10)
c
c Resultate umkopieren
c
        open(21,file='opt.res')
        open(22,file='opt.res.new')
 28     read(21,8000,END=29) line
        write(22,8000) line
        goto 28
 29     continue
c
c A/(A-1) Korrektur am tau
c
        tau=tau/2.
c
c Bestimmung der Reichweite des Potentials
c und der Schrittweite
c
	rp=sqrt(hb/(mt*omega))*3.
	dx=rp/npts
	print *,'dx = ',dx
c       
c Intergrationsroutine:
c Anfangsbedingungen:
c  Paritaet -
       v0=(pin**2-vopt(0.))
       psim  = (0.,0.)
       dpsim = sqrt(v0)
c  Paritaet +
       psip  = (1.,0.)
       dpsip = (0.,0.)
c Integration:
       do 100 i=1,npts
	  x=i*dx
c
c Runge Kutta Hilfsgroessen bestimmen
	  k1=dx*(vopt(x)-pin**2)*psip
	  k2=dx*(vopt(x+dx/2)-pin**2)*(psip+dx/2*dpsip+dx/8*k1)
	  k3=dx*(vopt(x+dx)-pin**2)*(psip+dx*dpsip+dx/2*k2)
c
c Interation
	  psip = psip+dx*(dpsip+(k1+2.*k2)/6.)
	 dpsip = dpsip+k1/6.+2./3.*k2+k3/6.
c
c Runge Kutta Hilfsgroessen bestimmen
	  k1=dx*(vopt(x)-pin**2)*psim
	  k2=dx*(vopt(x+dx/2)-pin**2)*(psim+dx/2*dpsim+dx/8*k1)
	  k3=dx*(vopt(x+dx)-pin**2)*(psim+dx*dpsim+dx/2*k2)
c
c Interation
	  psim = psim+dx*(dpsim+(k1+2.*k2)/6.)
	 dpsim = dpsim+k1/6.+2./3.*k2+k3/6.
100    continue
c
c Berechnung der logarithmischen Ableitungen:
c
       qm=dpsim/psim
       qp=dpsip/psip
       print *,qm,qp
c
c Berechnung von Reflektions+Transmissionsamplitude:
c
       pp = (pin*ic+qp)/(pin*ic-qp) * exp(-2.*ic*pin*rp) 
       pm = (pin*ic+qm)/(pin*ic-qm) * exp(-2.*ic*pin*rp) 
c
c Berechnung von t,r
c 
       r = (pp+pm)/2.
       t = (pp-pm)/2.
c
c noch A/(A-1) Korrektur
c
       r=r*2
       t=(t-1.)*2.+1.
        
       write(22,8007) tau,pin,t,r
       print *,'Reflektionsamplitude  : ',r,abs(r)**2
       print *,'Transmissionamplitude : ',t,abs(t)**2
       print *,'R+T = ',abs(r)**2+abs(t)**2 
c
c Ergebnissfiles umkopieren
c 
       close(21)
       close(22)
       open(21,file='opt.res.new')
       open(22,file='opt.res')
 2000  read(21,8000,END=2001) line
       write(22,8000) line
       goto 2000
 2001  continue

      
c
c Formats:
c
 8000  format(a100)
 8007  format(7e12.4)
       end


       complex function vopt(x)
       real x
C
C lokales:
       real pi
       real xe
C eingabedaten:
C mt:    Targetmasse
C mp:    Projetilmasse
C hb:    Planck-const
C omega: Oszillatorstaerke
C pin:   einlaufender Impuls
C tau:   elementare WW-Staerke
C nt:    Zahl der Targetzustaende
C np:    Zahl der Impulszustaende
       real mp,mt
       real hb
       real omega
       real pin
       complex tau
       integer nt,np
       
       common /param/mp,mt,hb,omega,pin,tau,nt,np ,pi
       
       if ((mt*omega/hb)*x**2.lt.16.) then       
	   xe=exp((-mt*omega/hb)*x**2)
       else
	   xe=0.
       endif
       vopt=2.*sqrt((mt*omega)/(pi*hb))*(2.*mp/hb)*tau*xe
       end




