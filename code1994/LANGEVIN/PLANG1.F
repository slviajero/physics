        program plang1
c------------------------------------------------------
c Dies ist ein experimentelles Programm zur Langevin
c Simulation im Phasenraum bei E>0.
c Es wird in einer Dimension gerechnet.
c --------------------------------------------------------
c Diese Version sampelt nicht die Anfangs- und Endpunkte
c durch den Gaussschen Zufallszahlengenerator, sondern
c bezieht Anfangs- und Endpunkt mit in die Langevin Simulation
c ein. Die Schleife ueber die Anfangskonfigurationen wird
c behalten, jedoch hat sie hier eine andere Bedeutung als
c in lang1 (evtl. Ergodizitaetsverletzung!)
c----------------------------------------------------------
c Testversion 1: teste woher Divergenzen des Messprozesses
c kommen um geeigneten Kern aufzubauen
c----------------------------------------------------------
c 
       implicit none
       integer nbmax
       parameter(nbmax=800)
       real db
       integer i,k1,j
       real pi
c
c externe Funktionen
c 
       real ran,rgau
       complex v,vr 
       external ran,rgau,v,vr
c
c Potentiale
c
       complex vp,vp1,vp2
       complex vk,vk1,vk2
       complex r,r2
c
c diverses
c
       real tt
       complex sv,sa,svv,sav
       complex ic
c
c Potentialanteil der Wirkung
c
       complex st,stv
c
c Langevin Schrittweiten
c
       real ds,dt
       integer nc,n0,nb
       integer ng
       integer l1,l2
c 
c Parameter des Potentials
c
       real m,b,v0,w,k,eta,v0k,a
c
c Anfangspunkte
c
       real xe,xa
       real pe,pa
c
c Anfangskonfigurationen
c 
       integer ia,na
c
c Ausgabeflags
c 
       logical cwrite
       logical swrite
       logical twrite

c
c Felder fuer Langevin 
c
c Pfad:
       complex x(0:nbmax+1), dx(0:nbmax+1)   
       complex p(0:nbmax+1), dp(0:nbmax+1)
c
c Noise:
       real dw(0:nbmax+1)
       complex z,zp
c
c Zufallsgenerator 
c
       double precision xr,xk
c
c globale variablen 
c
c      common/daten/x,dx,p,dp,dw
       COMMON/RG/Xr,XK
       common/const/pi
       common/par/m,b,v0,db,w,k,eta,v0k,a
c
c parameter
c
       pi     = 4.*atan(1.)
       ic     = (0.,1.)
       cwrite = .false.
       swrite = .true.
       twrite = .true.
c
c files
c
       open(10,file='test')
       open(12,file='error')
       if (swrite) open(15,file='act')
       open(16,file='conf')
       open(17,file='ratio')
c
c Hauptschleife
c 
1     continue
c
C  INITIALIZATION OF THE RANDOM NUMBER GENERATOR
      Xr=PI*1.D11
      XK=0.D00
      DO 10 K1=1,20
        TT=RAN()
10    CONTINUE
c
c Eingabe der Parameter 
c
      print *,'Modellparameter :'
      print *,'Zahl der beta-schritte :'
      read(*,*,END=9999) nb
      print *,'db :'
      read *,db
      print *,'Weite des Potentials :'
      read *,b
      print *,'Tiefe des Potentials :'
      read *,v0
c     print *,'Masse des Projectils :'
c     read *,m
      m=1.
c----------------------------------
c Vorsicht: m nicht implementiert!
c----------------------------------
      print *,'Energie :'
      read *,eta

      k=sqrt(2.*m*eta)
      w=2.*b*b
      a=2.*sqrt(2./pi)/b**(1./3.)
      v0k=v0*(2.*sqrt(2./pi)*b**2)**(4./3.)

      print *,'Langevin Parameter :'
      print *,'Zahl der gezogenen Konfigurationen :'
      read *, nc
      print *,'Gleichgewichtskonfigurationen      :'
      read *,ng
      print *,'Zahl der ausgelassenen Schritte    :'
      read *,n0
      print *,'Laenge des Zeitintervalls :'
      read *,dt
      print *,'Simulationszeit :'
      print *,nc*n0*dt
      
c      print *,'Zahl der Anfangskonfigurationen :'
c      read *,na
      na=1 

      write(15,*) na,nb,nc,dt
      write(15,*) m,b,v0,db,eta
      ds=sqrt(dt*2.)

      sa=0.
      sav=0.
      do 1000 ia=1,na           
      sv=0.
      svv=0.
c
c Anfangspunkte
c
c     print *,'Anfangswerte samplen :'
      do 20 j=1,3
             xe=rgau()*b
20    continue
      do 21 j=1,3
             xa=rgau()*b
21    continue
      pa=0.
      pe=0.
c
c Anfangswert, straight line zwischen Anfangs und Endpunkt
c
c       print *,'Anfangswert fuer Langevin Gleichung machen :'
        do 100 i=0,(nb+1)
              x(i)=cmplx((xe-xa)*i/real(nb+1)+xa)
	      p(i)=cmplx(0.)
100     continue
c
c Langevin Simulation
c 
c
c Gleichgewichtslauf
c
c         print *,'Gleichgewichtslauf :',ia
          do 190 l2=1,ng 
c          print *,l2
c
c ein schritt
c
             call makran(dw,nb,ds)
             do 191 i=0,nb+1
                call s(i,x,p,z,zp,nb)
                dx(i)=z*dt+dw(i)
		dp(i)=zp*dt
191          continue
             do 193 i=0,nb+1
                x(i)=x(i)+dx(i)
		p(i)=p(i)+dp(i)
193          continue   
c
c ende des schrittes
c
190     continue
c
c eigentliche simulation
c
199     continue
c       print *,'Simulationsbeginn ',ia
        do 200 l1=1,nc
c          if (mod(l1,1000).eq.0) print *,l1
c
c auslassung (n0 schritte)
c
          do 210 l2=1,n0 
c
c ein schritt
c 
             call makran(dw,nb,ds)
             do 220 i=0,nb+1
                call s(i,x,p,z,zp,nb)
                dx(i)=z*dt+dw(i)
		dp(i)=zp*dt
220          continue
             do 240 i=0,nb+1
                x(i)=x(i)+dx(i)
		p(i)=p(i)+dp(i)
240          continue   
c
c ende des schrittes
c   

        if ((mod(l1,10).eq.0)) then
             write(10,*) dt*l1, real(x(nb/2)),aimag(x(nb/2)),
     A                          real(p(nb/2)),aimag(p(nb/2))
	endif
210     continue
c
c Potentialanteil der Wirkung berechnen
c
        st  = cmplx(0.)
	stv = cmplx(0.)
        do i=0,nb
	  r   = x(i)
          r2  = r**2
          vp  = v0*exp(-r2/w)
          vp1 =-r/b**2*vp
          vp2 =-1/b**2*(1.-r2/b**2)*vp
	  vk  = v0k/(r2+a**2)**2
	  vk1 = v0k*(-4.)*r/(r2+a**2)**3
	  vk2 = v0k*(4.)*(7.*r2-a**2)/(r2+a**2)**4
          st  = st 
     A        + 2.*(.5*p(i)**2 - eta)*(vp - vk)
     B        + vp**2 - vk**2
     C        + ic*p(i)*(vp1 - vk1)
     D        - .5*(vp2 - vk2)
        end do
        stv=exp(st*db)
c
c Korrektur der Potentialform
c
c Anfangspunkt:
c
        r  = x(0)
	r2 = r**2
   	vp  = v0*exp(-r2/w)
        vk  = v0k/(r2+a**2)**2
	stv=stv*(v0k/v0)
c
c Endpunkt:
c
        r  = x(nb+1)
	r2 = r**2
   	vp  = v0*exp(-r2/w)
        vk  = v0k/(r2+a**2)**2
	stv=stv*(v0k/v0)
c
c und Rausschreiben:
c
        sv=sv+stv
        svv=svv+stv**2
        if (twrite.and.(mod(l1,100).eq.0)) then
            write(15,*) l1*dt,real(st),aimag(st),
     A                   real(exp(st)),aimag(exp(st))
        endif

c
c Konfiguration auf file schreiben
c
200    continue
c
c hier endet die Simulationsschleife ueber nc Konfigurationen
c
c berechne mittleren potentialanteil der Wirkung
c fuer eine anfangskonfiguration
       sv=sv/nc
       svv=(svv/nc-sv**2)
       sa=sa+sv
       sav=sav+svv
1000   continue
c
c hier endet die Simulationschleife ueber na Anfangs
c konfigurationen
       sa=sa/na
       sav=sav/na
c
c letzte Konfiguration abspeichern
c
	if (cwrite) then
          close(16)
	  open(16,file='conf')
	  do i=0,nb+1
	     write(16,1010) i,x(i)
          end do
	  write(16,*) xr,xk
        end if
c
c  Verhaeltniss ausgeben
c
       write(17,*) (nb+1)*db, 1/((nb+1)*db), sa
c      write(12,*) (nb+1)*db, 1/((nb+1)*db), sav, sqrt(sav)
       goto 1
1010   format(i3,e12.4)
9999   continue
       end
      
      subroutine s(i,x,p,z,zp,nb)
c
c berechnet die ableitung der Wirkung am iten Schritt
c fuer ein gausspotential der form
c v0*exp(-x**2 /(2*b**2))
c
      implicit none
      integer nbmax
      parameter(nbmax=800)
      real b,v0,db,m,pi,w,k,eta,v0k,a
      integer i,nb
      complex v
      external v
      complex x(0:nbmax+1)
      complex p(0:nbmax+1)
      complex z,zp
      complex ic
      complex vp1,vp2,vp3,r,r2,vp
      common/const/pi
      common/par/m,b,v0,db,w,k,eta,v0k,a

      ic=(0.,1.)

      r=x(i)
      r2=r**2
c
c Ableitungen des Potentials
c
      vp  = v0*exp(-r2/w)
      vp1 =-v0*r/b**2*exp(-r2/w)
      vp2 =-v0/b**2*(1.-r2/b**2)*exp(-r2/w)
      vp3 = v0*r/b**4*(3.-r2/b**2)*exp(-r2/w)

      if (i.eq.0) then
      
          z=-(  (2.*(.5*p(i)**2 + vp - eta)*vp1 
     A           + ic*p(i)*vp2 
     B           - .5*vp3 
     C          )*db
     D         + ic*(p(i)-k)
     E         + r/b**2        
     F        )
          zp=-(  (2.*(.5*p(i)**2 + vp - eta)*p(i)
     A            + ic*vp1
     B           )*db
     C         -ic*(x(i+1)-x(i))
     D        )
      
      elseif (i.eq.(nb+1)) then
      
          z=-(  
     A         - ic*(p(i-1)-k)
     B         + r/b**2
     D       )
          zp= 0. 
      
      else

          z=-(  (2.*(.5*p(i)**2 + vp - eta)*vp1 
     A           + ic*p(i)*vp2 
     B           - .5*vp3 
     C          )*db
     D         + ic*(p(i)-p(i-1)) 
     E       )
          zp=-(  (2.*(.5*p(i)**2 + vp - eta)*p(i)
     A            + ic*vp1
     B           )*db
     C         -ic*(x(i+1)-x(i))
     D        )

      endif
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


      real function rgau()
      real x1,x2
      real ran
      external ran
      real pi
      common/const/pi
      x1=ran()
      x2=ran()
      rgau=sqrt(-2.*alog(x1))*cos(2.*pi*x2)
      end


      subroutine makran(w,n,ds)
      integer nbmax
      parameter(nbmax=800)
      real w(0:nbmax+1)
      real x1,x2
      real cs,sn,lg
      real ds
      integer n
      integer i
      real pi
      common/const/pi
c
c macht die zufallszahlen fuer einen schritt
c
      do 100 i=0,n,2
10         x1=ran() 
	   if (x1.lt.1e-16) goto 10
           x2=ran()
           cs=cos(2.*pi*x2)
           sn=sin(2.*pi*x2)
           lg=sqrt(-2.*alog(x1))*ds
           w(i)  =lg*cs
           w(i+1)=lg*sn
100   continue

      end

