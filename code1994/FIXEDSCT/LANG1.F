       program langv1
c--------------------------------------------------------
c Fuehrt Langevin Simulation durch um Konfigurationen zu
c gewinnen, die nach der Wirkung von einem Teilchen im
c Gausspotential verteilt sind. Berechnet wird direkt die 
c Streulaenge.
c --------------------------------------------------------
       integer ncmax,nbmax
       parameter(ncmax=100,nbmax=100)
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
       real tt,st,sv,sa
       real ds,dt
       integer nl,nc,n0,nb
       integer l1,l2
       real m,b,v0,w
       real xi(3),xf(3)
       integer ia,na
       real r
       logical cwrite
c
c Felder fuer Langevin 
c
       real x(3,0:nbmax+1), dx(3,0:nbmax+1)       
       real dw(3,nbmax+1)
       real z(3)
c
c Zufallsgenerator 
c
       double precision xr,xk
c
c globale variablen 
c
       COMMON/RG/Xr,XK
       common/const/pi
       common/par/m,b,v0,db,w
c
c parameter
c
       pi =4.*atan(1.)
       cwrite=.false.
c
c files
c
       open(10,file='test')
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
      read *,nb
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
      print *,'Zahl der ausgelassenen Schritte    :'
      read *,n0
      print *,'Laenge des Zeitintervalls :'
      read *,dt
      print *,'Simulationszeit :'
      print *,nc*n0*dt
      print *,'Korrelationszeit :'
      print *,n0*dt
      
      print *,'Zahl der Anfangskonfigurationen :'
      read *,na
       
      write(10,*) na,nb,nc,dt
      write(10,*) m,b,v0,db
      ds=sqrt(dt*2.)

      sa=0.
      do 1000 ia=1,na           
      sv=0.
c
c Anfangspunkte
c
      print *,'Anfangswerte samplen :'
      do 20 j=1,3
             xi(j)=rgau()*b
20    continue
      do 21 j=1,3
             xf(j)=rgau()*b
21    continue
c
c Anfangswert, straight line zwischen Anfangs und Endpunkt
c
        print *,'Anfangswert fuer Langevin Gleichung machen :'
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
          print *,'Gleichgewichtslauf :',ia
          do 190 l2=1,n0 
c          print *,l2
c
c ein schritt
c
             call makran(dw,nb,ds)
             do 191 i=1,nb
                call s(i,x,z)
                do 192 j=1,3
                  dx(j,i)=z(j)*dt+dw(j,i)
c                  print *,i,j,z(j)*dt
192             continue
191          continue
             do 193 i=1,nb
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
195          continue
           
190     continue
c
c eigentliche simulation
c
        print *,'Simulationsbeginn ',ia
        do 200 l1=1,nc
c
c auslassung (n0 schritte)
c
          do 210 l2=1,n0 
c
c ein schritt
c 
             call makran(dw,nb,ds)
             do 220 i=1,nb
                call s(i,x,z)
                do 230 j=1,3
                  dx(j,i)=z(j)*dt+dw(j,i)
230             continue
220          continue
             do 240 i=1,nb
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
         do 400 i=0,nb
            r=0.
            do 410 j=1,3
              r=r+x(j,i)*x(j,i)
410         continue
            r=sqrt(r) 
            st=st+v(r)*db
400     continue
        sv=sv+exp(st)           
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
c berechne mittleren potentialanteil der Wirkung
c fuer eine anfangskonfiguration
       sv=sv/nc
       sa=sa+g0(xf,xi,db*nb)/sv       
1000   continue
c
c Normierung korrigieren
c
c 
       sa=sa/na
       sa=sa*v0**2*b**6*(2.*pi)**3
       sa=sa*(2.*m*(2.*pi)**3/(4.*pi))**2/(2.*pi)**3       
       sa=sa*(nb*db/(2.*pi*m))**(3./2.)
       
       print *,'Steulaenge  : ',sqrt(sa)
       print *,'Querschnitt : ',4*pi*sa

1010   format(3e12.4)
       end
      
      subroutine s(i,x,z)
c
c berechnet die ableitung der Wirkung am iten Schritt
c
      integer nbmax
      parameter(nbmax=100)
      real b,v0,db,m,pi,w
      integer i,j
      real x(3,0:nbmax+1)
      real z(3)
      real vp,u,r
      external vp
      common/const/pi
      common/par/m,b,v0,db,w

      r=0.
      do 10 j=1,3
         r=r+x(j,i)*x(j,i)
10    continue
      r=sqrt(r)     
      u=vp(r)*db
      do 100 j=1,3
        z(j)=m/db*(x(j,i+1)-2.*x(j,i)+x(j,i-1))-u
100   continue
      end


      real function vp(r)
c
c ableitung des potentials (radialsymetrisch)
c
      real pi,m,b,v0,db,w 
      real r 
      common/const/pi
      common/par/m,b,v0,db,w

      vp=-2.*r/w*v0*exp(-r*r/w)
      end


      real function v(r)
c
c potential (radialsymetrisch)
c
      real pi,m,b,v0,db,w 
      real r 
      common/const/pi
      common/par/m,b,v0,db,w

      v=v0*exp(-r*r/w)
      end


      real function g0(x,y,beta)
c
c berechnet freie euklidische Greensfunktion
c
      real pi,m,b,v0,db,w 
      integer j
      real x(3),y(3)
      real r
      common/const/pi
      common/par/m,b,v0,db,w      
      
      r=0.
      do 10 j=1,3
         r=r+(x(j)-y(j))**2
10    continue 
             
      g0=(m/(2.*pi*beta))**(3./2.)*exp(-m/(2.*beta)*r)

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
      parameter(nbmax=100)
      real w(3,1:nbmax+1)
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
      do 100 i=1,n,2
         do 100 j=1,3
           x1=ran()
           x2=ran()
           cs=cos(2.*pi*x2)
           sn=sin(2.*pi*x2)
           lg=sqrt(-2.*alog(x1))*ds
           w(j,i)  =lg*cs
           w(j,i+1)=lg*sn
100   continue

      end
