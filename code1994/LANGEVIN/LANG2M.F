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
       integer ncmax,nbmax
       parameter(ncmax=100,nbmax=800)
       real db
       integer n,i,k1,j
       real pi
c
c externe Funktionen
c 
       real ran,rgau,v,vp
       external ran,rgau,v,vp
c
c
c
       real tt,st,sv,sa,sb,svv
       real samin,samax
       real ds,dt
       integer nl,nc,n0,nb
       integer ng
       integer l1,l2
       real m,b,v0,w
       real xi(3),xf(3)
       integer ia,na
       real r
       logical cwrite
       logical swrite
c
c Felder fuer Langevin 
c
       real x(3,0:nbmax+1), dx(3,0:nbmax+1)       
       real dw(3,0:nbmax+1)
       real r2(0:nbmax+1)
       real z(3)
c
c Zufallsgenerator 
c
       double precision xr,xk
c
c globale variablen 
c
c      common/daten/x,dx,dw,r2
       COMMON/RG/Xr,XK
       common/const/pi
       common/par/m,b,v0,db,w
c
c parameter
c
       pi =4.*atan(1.)
       cwrite=.false.
       swrite=.false.
c
c files
c
       open(10,file='test')
       open(11,file='res')
       open(12,file='error')
       open(14,file='relax')
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
      print *,'Masse des Projectils :'
      read *,m

      w=2.*b*b

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

      write(10,*) na,nb,nc,dt
      write(10,*) m,b,v0,db
      ds=sqrt(dt*2.)

      sa=0.
      do 1000 ia=1,na           
      sv=0.
      svv=0.
c
c Anfangspunkte
c
c     print *,'Aufsetzen auf Konfiguration ?'
c     read *,i
      i=0
      if (i.eq.1) then 
	 rewind(16)
	 do i=0,nb+1
	    read(16,1010) n,(x(j,i),j=1,3)
         end do
	 read(16,*) xr,xk
	 goto 199
      endif
c     print *,'Anfangswerte samplen :'
      do 20 j=1,3
             xi(j)=rgau()*b
20    continue
      do 21 j=1,3
             xf(j)=rgau()*b
21    continue
c
c Anfangswert, straight line zwischen Anfangs und Endpunkt
c
c       print *,'Anfangswert fuer Langevin Gleichung machen :'
        do 100 i=0,(nb+1)
           do 110 j=1,3
              x(j,i)=(xf(j)-xi(j))*i/real(nb+1)+xi(j)
110        continue
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
                call s(i,x,z,nb,r2)
                do 192 j=1,3
                   dx(j,i)=z(j)*dt+dw(j,i)
c                  print *,i,j,z(j)*dt
192             continue
191          continue
             do 193 i=0,nb+1
                do 194 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
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
                call s(i,x,z,nb,r2)
                do 230 j=1,3
                  dx(j,i)=z(j)*dt+dw(j,i)
230             continue
220          continue
             do 240 i=0,nb+1
                do 250 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
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
            r=r2(i) 
            st=st+v(r)*db
400     continue 
	st=st+v(r2(nb+1))*db/2
	st=st+v(r2(0))/2*db

        if (swrite) write(15,*) l1,exp(st)
        sv=sv+exp(st) 
	svv=svv+exp(2*st)
c
c Konfiguration auf file schreiben
c
        if (.false.) then
	  write(10,*) 'Konfiguration :',l1
          do 300 i=0,nb+1
              write(10,1010) i,(x(j,i),j=1,3)
	      write(10,1010) i,sqrt(r2(i))
	      write(10,1010) i,(dw(j,i),j=1,3)
	      write(10,1010) i,(dx(j,i),j=1,3)
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
1000   continue
c
c hier endet die Simulationschleife ueber na Anfangs
c konfigurationen
       sa=sa/na
c
c letzte Konfiguration abspeichern
c
	close(16)
	open(16,file='conf')
	do i=0,nb+1
	   write(16,1010) i,(x(j,i),j=1,3)
        end do
	write(16,*) xr,xk
c
c  Verhaeltniss ausgeben
c
       write(17,*) (nb+1)*db, sa

c
c hier wird der Bornanteil dazugegeben
c
       samin=gborn((nb+1)*db)/(sa+sqrt(svv/nc))
       samax=gborn((nb+1)*db)/(sa-sqrt(svv/nc))
       sa   =gborn((nb+1)*db)/(sa)

c
c Normierung korrigieren
c
c
       sb=gborn((nb+1)*db)
       sb=sb*(2.*m*(2.*pi)**3/(4.*pi))**2/(2.*pi)**3       
       sb=sb*((nb+1)*db/(2.*pi*m))**(3./2.)
       sa=sa*(2.*m*(2.*pi)**3/(4.*pi))**2/(2.*pi)**3       
       sa=sa*((nb+1)*db/(2.*pi*m))**(3./2.)
      

       samax=samax*(2.*m*(2.*pi)**3/(4.*pi))**2/(2.*pi)**3       
       samax=samax*((nb+1)*db/(2.*pi*m))**(3./2.)


       samin=samin*(2.*m*(2.*pi)**3/(4.*pi))**2/(2.*pi)**3       
       samin=samin*((nb+1)*db/(2.*pi*m))**(3./2.)


       print *,'Steulaenge  : ',sqrt(sa)
       print *,'Querschnitt : ',4*pi*sa
       print *,'Born        : ',4*pi*sb

          write(11,*) (nb+1)*db,1/((nb+1)*db),4*pi*sa,4*pi*sb
	  write(12,*) (nb+1)*db,4*pi*sa,4*pi*samin,4*pi*samax
          write(14,*) ng*dt,nc*dt,4*pi*sa,4*pi*sb


c         write(11,*) ng*dt , 4*pi*sa
       goto 1
1010   format(i3,3e12.4)
9999   continue
       end
      
      subroutine s(i,x,z,nb,r2)
c
c berechnet die ableitung der Wirkung am iten Schritt
c
      integer nbmax
      parameter(nbmax=800)
      real b,v0,db,m,pi,w
      integer i,j,nb
      real x(3,0:nbmax+1)
      real r2(0:nbmax+1)
      real z(3)
      real vp,u,r
      common/const/pi
      common/par/m,b,v0,db,w
 
      r=0.
      do 10 j=1,3
         r=r+x(j,i)*x(j,i)
10    continue
      r2(i)=r
      if (r/w.lt.40) then
          vp =-2./w*v0*exp(-r/w)*db
      else
	  vp = 0.
      endif
      if (i.eq.0) then
        do 100 j=1,3
          z(j)=m/db*(x(j,i+1)-x(j,i))
     A          -vp*x(j,i)/2.-x(j,i)/b**2
100     continue
      elseif (i.eq.(nb+1)) then
        do 110 j=1,3
          z(j)=m/db*(-x(j,i)+x(j,i-1))
     A          -vp*x(j,i)/2.-x(j,i)/b**2
110     continue
      else
        do 120 j=1,3
          z(j)=m/db*(x(j,i+1)-2.*x(j,i)+x(j,i-1))-vp*x(j,i)
120     continue
      endif
      end


      real function v(r)
c
c potential (radialsymetrisch)
c r: quadrat von r
c
      real pi,m,b,v0,db,w 
      real r 
      common/const/pi
      common/par/m,b,v0,db,w

      v=v0*exp(-r/w)
      end


      real function gborn(beta)
c
c berechnet erste Ordnung Bornapproximation
c
      real beta
      real m,b,v0,db,w
      real pi
      common/const/pi
      common/par/m,b,v0,db,w      
      
             
      gborn=v0**2*b**3/sqrt(1+beta/m/w)**3*pi**(1.5)

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
      real w(3,0:nbmax+1)
      real x1,x2,g1,g2
      real cs,sn,lg
      real ds
      integer n
      integer i,j
      real pi
      common/const/pi
c
c macht die zufallszahlen fuer einen schritt
c
      do 100 i=0,n,2
         do 100 j=1,3
10         x1=ran() 
	   if (x1.lt.1e-16) goto 10
           x2=ran()
           cs=cos(2.*pi*x2)
           sn=sin(2.*pi*x2)
           lg=sqrt(-2.*alog(x1))*ds
           w(j,i)  =lg*cs
           w(j,i+1)=lg*sn
100   continue

      end
