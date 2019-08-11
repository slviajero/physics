       program langv2           
c--------------------------------------------------------
c Fuehrt Langevin Simulation durch um Konfigurationen zu
c gewinnen, die nach der Wirkung von einem Teilchen im
c Gausspotential verteilt sind. Berechnet wird direkt die 
c Streulaenge.
c --------------------------------------------------------
c Diese Version sampelt nicht die Anfangs- und Endpunkte
c durch den Gaussschen Zufallszahlengenerator, sondern
c bezieht Anfangs- und Endpunkt mit in die Langevin Simulation
c ein. Die Schleife ueber die Anfangskonfigurationen wird
c behalten, jedoch hat sie hier eine andere Bedeutung als
c in lang1
c----------------------------------------------------------
c oszillator version des langevin programms
c----------------------------------------------------------
c diese version implementiert die greensfunktion des
c oszillators exakt, dh. eine wirkung mit cosh und sinh
c termen wird verwendet
c----------------------------------------------------------
c zusaetzlich wird mit gefalteten potential verglichen
c ausgegeben wird nurmehr das verhaeltnis der querschnitte 
c-----------------------------------------------------------
       implicit none 
       integer ncmax,nbmax
       parameter(ncmax=100,nbmax=800)
       real db
       integer i,k1,j
       real pi
c
c externe Funktionen
c 
       real ran,rgau,v,vf
       external ran,rgau,v,vf
c
c
c
       real tt,st,sv,sa,svv,sav
       real ds,dt
       integer nc,n0,nb
       integer ng
       integer l1,l2
c Parameter des Projectils:
       real m,b,v0,w,mb
c Parameter des Oszi:
       real omega,mo,xo,mob,f1,f2
c Anfangskonfiguration:
       real xi(3),xf(3)
       real qi(3),qf(3)
c zahl der Anfangskonfigs:
       integer ia,na
       logical cwrite
       logical twrite
c
c Felder fuer Langevin 
c      
       real r2(0:nbmax+1), r2t(0:nbmax+1) 
       real x(3,0:nbmax+1), dx(3,0:nbmax+1)
       real q(3,0:nbmax+1), dq(3,0:nbmax+1)       
       real dw(3,0:nbmax+1)
       real dwq(3,0:nbmax+1)
       real z(3)
       real zq(3)
c
c Zufallsgenerator 
c
       double precision xr,xk
c
c globale variablen 
c
       COMMON/RG/Xr,XK
       common/const/pi
       common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2
c
c parameter
c
       pi =4.*atan(1.)
       cwrite=.false.
       twrite=.true.
c
c files
c
       open(10,file='test')
       open(12,file='error')
       open(14,file='ratio')

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
      print *,'Masse des Projectils :'
      read *,m


      print *,'Oszillatorstaerke :'
      read *,omega
      print *,'Oszi Masse:'
      read *,mo


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
      print *,'Korrelationszeit :'
      print *,n0*dt
      
c      print *,'Zahl der Anfangskonfigurationen :'
c      read *,na
      na=1 
      
      write(10,*) 'Lang2hof: Version 1.12.1993'
      write(10,*) na,nb,nc,dt
      write(10,*) m,b,v0,db
      write(10,*) mo,omega
c
c berechne bestimmte kombinationen von werten im voraus
c
      ds=sqrt(dt*2.)
      xo=sqrt(mo*omega)
      w=2.*b*b
      mob=mo/db
      mb=m/db
      f1=mo*omega/sinh(db*omega)
      f2=cosh(db*omega)


      sa=0.
      sav=0.
      do 1000 ia=1,na           
      sv=0.
      svv=0.
c-----------------------------------------------------------
c Anfangspunkte
c
c     print *,'Anfangswerte samplen :'
      do 23 j=1,3
             qi(j)=rgau()*xo
 23   continue
      do 24 j=1,3
             qf(j)=rgau()*xo
 24   continue
      do 20 j=1,3
             xi(j)=rgau()*b+qi(j)
20    continue
      do 21 j=1,3
             xf(j)=rgau()*b+qf(j)
21    continue
c--------------------------------------------------------------
c Anfangswert, straight line zwischen Anfangs und Endpunkt
c das kann man sicher besser machen !
c
c       print *,'Anfangswert fuer Langevin Gleichung machen :'
        do 100 i=0,(nb+1)
           do 110 j=1,3
              x(j,i)=(xf(j)-xi(j))*i/real(nb+1)+xi(j)
              q(j,i)=(qf(j)-qi(j))*i/real(nb+1)+qi(j)
110        continue
100     continue
c-------------------------------------------------------------
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
             call makran(dw,dwq,nb,ds)
             do 191 i=0,nb+1
                call s(i,x,q,z,zq,nb,r2,r2t)
                do 192 j=1,3
                   dx(j,i)=z(j)*dt+dw(j,i)
                   dq(j,i)=zq(j)*dt+dwq(j,i)
c                  print *,i,j,z(j)*dt
192             continue
191          continue
             do 193 i=0,nb+1
                do 194 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
                  q(j,i)=q(j,i)+dq(j,i)
194             continue
193          continue   
c
c ende des schrittes
c
c             do 195 i=1,nb
c                do 195 j=1,3
c                   print *,i,j,x(j,i),dx(j,i),dw(j,i)
c195          continue
           
190     continue
c
c eigentliche simulation
c
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
             call makran(dw,dwq,nb,ds)
             do 220 i=0,nb+1
                call s(i,x,q,z,zq,nb,r2,r2t)
                do 230 j=1,3
                  dx(j,i)=z(j)*dt+dw(j,i)
                  dq(j,i)=zq(j)*dt+dwq(j,i)
230             continue
220          continue
             do 240 i=0,nb+1
                do 250 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
                  q(j,i)=q(j,i)+dq(j,i)
250             continue
240          continue   
c
c ende des schrittes
c
210     continue
c
c Potentialanteil der Wirkung berechnen
c
         st=0.
         do 400 i=1,nb
            st=st+(v(r2(i))-vf(r2t(i)))*db
400     continue
	st=st+(v(r2(nb+1))-vf(r2t(nb+1)))*db/2
	st=st+(v(r2(0))-vf(r2t(0)))*db/2

c        print *,'Wirkung der ',l1,'.  Konf.:',st
c        write(10,*) l1,exp(st)
        sv=sv+exp(st) 
	svv=svv+exp(2*st)
	if ((mod(l1,100).eq.0).and.twrite) then
	   write(10,*) dt*l1,sv/l1
        endif

c
c Konfiguration auf file schreiben
c
        if (cwrite) then
          do 300 i=0,nb+1
              write(10,1010) (x(j,i),j=1,3)
300       continue
        endif
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
       write(14,*) (nb+1)*db,1/((nb+1)*db),sa
       write(12,*) (nb+1)*db,1/((nb+1)*db),sav,sqrt(sav)
       goto 1
1010   format(3e12.4)
9999   continue
       end
      
      subroutine s(i,x,q,z,zq,nb,r2,r2t)
c
c berechnet die ableitung der Wirkung am iten Schritt
c
      integer nbmax
      parameter(nbmax=800)
      real b,v0,db,m,pi,w,mb,mob
      integer i,j,nb
      real r2(0:nbmax+1),r2t(0:nbmax+1)
      real x(3,0:nbmax+1),q(3,0:nbmax+1)
      real z(3),zq(3)
      real vp,r,rt
      real omega,mo,xo
      common/const/pi
      common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2
 
      r=0. 
      rt=0.
      do 10 j=1,3
         r=r+(x(j,i)-q(j,i))**2
	 rt=rt+x(j,i)**2
10    continue
      r2(i)=r
      r2t(i)=rt
      vp =-2./w*v0*exp(-r/w)*db
      if (i.eq.0) then
        do 100 j=1,3
c Projectile: kinetic
          z(j)=mb*(x(j,i+1)-x(j,i))
c             potential
     A          -vp*(x(j,i)-q(j,i))/2.
c             endpoint
     B          -(x(j,i)-q(j,i))/b**2
c Target: kinetic
          zq(j)=f1*(q(j,i+1)-f2*q(j,i))
c         potential
     B         -vp*(q(j,i)-x(j,i))/2.
c         endpoint
     C         -xo**2*q(j,i)-(q(j,i)-x(j,i))/b**2
100     continue
      elseif (i.eq.(nb+1)) then
        do 110 j=1,3
c Projectile: kinetic
          z(j)=mb*(-x(j,i)+x(j,i-1))
c             potential
     A          -vp*(x(j,i)-q(j,i))/2.
c             endpoint
     B          -(x(j,i)-q(j,i))/b**2
c Target: kinetic
          zq(j)=f1*(-f2*q(j,i)+q(j,i-1))
c         potential
     B         -vp*(q(j,i)-x(j,i))/2.
c         endpoint
     C         -xo**2*q(j,i)-(q(j,i)-x(j,i))/b**2
110     continue
      else
        do 120 j=1,3
c Projectile: kinetic    
          z(j)=mb*(x(j,i+1)-2.*x(j,i)+x(j,i-1))
c             potential
     A         -vp*(x(j,i)-q(j,i))
c Target: kinetic
          zq(j)=f1*(q(j,i+1)-2.*f2*q(j,i)+q(j,i-1))
c         potential
     B         -vp*(q(j,i)-x(j,i))
120     continue
      endif
      end


      real function v(r)
c
c potential (radialsymetrisch)
c r ist der radius quadrat !!!
c
      implicit none
      real pi,m,b,v0,db,w
      real omega,mo,xo,mb,mob,f1,f2
      real r 
      common/const/pi
      common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2

      v=v0*exp(-r/w)
      end

      
      real function vf(r)
c
c potential (radialsymetrisch) gefaltet mit grundzustand
c des harmonischen oszi
c r ist der radius quadrat !!!
c
      implicit none
      real pi,m,b,v0,db,w
      real omega,mo,xo,mb,mob,f1,f2
      real r,alpha 
      common/const/pi
      common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2

      alpha=w*mo*omega
      vf=v0*sqrt(alpha/(1+alpha))**3*exp(-r/w*alpha/(1+alpha))

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


      subroutine makran(w,wq,n,ds)
      implicit none
c
c maximale zahl der beta schritte
c
      integer nbmax
      parameter(nbmax=800)
c
c felder fuer zufallszahlen
c
      real w(3,0:nbmax+1)
      real wq(3,0:nbmax+1)
      real x1,x2
      real cs,sn,lg
      real ds
      real ran
      external ran
      integer n
      integer i,j
      real pi
      common/const/pi
c
c macht die zufallszahlen fuer einen schritt
c fuer das projektil
c
      do 100 i=0,n,2
         do 100 j=1,3
101        x1=ran()
	   if (x1.lt.1e-16) goto 101
           x2=ran()
           cs=cos(2.*pi*x2)
           sn=sin(2.*pi*x2)
           lg=sqrt(-2.*alog(x1))*ds
           w(j,i)  =lg*cs
           w(j,i+1)=lg*sn
100   continue
c
c zufallszahlen fuer harmonischen oszillator
c
      do 110 i=0,n,2
         do 110 j=1,3
111        x1=ran()
	   if (x1.lt.1e-16) goto 111
           x2=ran()
           cs=cos(2.*pi*x2)
           sn=sin(2.*pi*x2)
           lg=sqrt(-2.*alog(x1))*ds
           wq(j,i)  =lg*cs
           wq(j,i+1)=lg*sn
110   continue
      end
