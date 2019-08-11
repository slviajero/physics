       program optx
       implicit none
c parameter:
       integer npts,np,nbeta
       real*8 pmax,pmin
       real*8 betama,betami
       parameter(npts=8000,np=400,nbeta=80)
       parameter(pmax=4.,pmin=0.)
       parameter(betama=40.0,betami=0.)
C
C npts: Zahl der Integrationspunkte
C-----------------------------------------------------------
C Autor:   Stefan Lenz
C Datum:
C Version: 0.0
C-----------------------------------------------------------
C Berechnet < \rho V V > fuer beliebige Potential
C Diese Version benutzt die Runge Kutta Methode aus Abramowitz
C S897
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
       complex*16 qp,psinor
       real*8 p
       real*8 beta
       real*8 ephase
       real*8 resp(0:np),resb(nbeta) 
       complex*16 sum
       real *8 sumb
       integer ip,ib,nb
       real*8 db,sa,dp
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

       common /param/m,b,v0,pi,p

       pi = 4.D00*atan(1.D00)
       ic = (0. , 1.D00)
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
c Files oeffnen
c
       open(10,file='integ')
       open(11,file='beta')

c
c Bestimmung der Reichweite des Potentials
c und der Schrittweite
c
	rp=b*10.
	dx=rp/npts
	dp=(pmax-pmin)/np
	do 200 ip=0,np 
	p=dp*ip
c       
c Intergrationsroutine:
c Anfangsbedingungen:
       psip  = (0.,0.)
       dpsip = sqrt(-vopt(0.D00)) 
c Integration:
       sum=(0.,0.)
       do 100 i=1,npts
	  x=i*dx
	  sum=sum+psip*(vopt(x)+p**2)*dx*x
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
c Berechnung der Streuphase:
c
       ephase=-exp(-2*ic*p*rp)*(ic*p+qp)/(ic*p-qp)
c
c Berechnung der Normierung
c
       psinor= ( -exp(-ic*p*rp)/(ic*p-qp))/psip
c
c Ummormieren der Intergrals
c
       sum=sum*psinor/(2*pi)**(3./2.)/(2*m)
       resp(ip)=abs(sum)**2
200    continue

c
c Laplace Trafo
c
       db=(betama)/nbeta
       do ib=1,nbeta
	  sumb=(0.,0.)
	  beta=db*ib
	  do ip=0,np
	     p=dp*ip
        sumb=sumb+resp(ip)*exp(-beta*p**2/(2.D00*m))*(p**2)*dp
          end do
          resb(ib)=(4*pi)**3*sumb   
       end do
c
c Ausgabe (mit richtigem Faktor)  
c
        do nb=1,nbeta
	sa=resb(nb)
        sa=sa*(2.*m*(2.*pi)**3/(4.*pi))**2/(2.*pi)**3       
        sa=sa*(nb*db/(2.*pi*m))**(3./2.)
        
        print *,'Steulaenge  : ',sqrt(sa)
        print *,'Querschnitt : ',4*pi*sa

        write(11,*) nb*db, 1/(nb*db) , 4*pi*sa
        end do
 

       do ip=1,np
	 sa=(4*pi)**2*resp(ip)*(2.*pi)**3
         write(10,*) ip*dp,sa
       end do

        end




       complex*16 function vopt(x)
       implicit none
       real*8 x
C
C lokales:
       real*8 pi
       real*8 m,b,v0,p,xe
       
       common /param/m,b,v0,pi,p
       
       if ((x/b)**2.lt.80.) then       
           xe=exp(-(x/b)**2/2.)
       else
           xe=0.D00
       endif
c      if (x.lt.b) then
c xe=1.0
c      else
c xe=0.
c      endif
       vopt=2.*v0*m*xe-p**2
       end

