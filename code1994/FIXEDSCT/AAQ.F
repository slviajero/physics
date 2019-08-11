       program ntest
       integer ncmax,nbmax,namax
       parameter(namax=100,ncmax=1000,nbmax=10)
       real db
       integer n,i,k1,j
       real pi
       real tt,tm,hm,sm,ids
       real v
       external v
       real r
       real ds,dt
       integer nl,nc,n0,nb
       integer l1,l2
       real m,b,v0,w
       real xi(3),xf(3)
       integer ia,na
c
c Felder 
c
       real x(3,0:nbmax+1,ncmax), mu(3,nbmax), sig(3,nbmax)
       real z(3)
c
c parameter 
c
       double precision xr,xk
       COMMON/RG/Xr,XK
       common/const/pi
       common/par/m,b,v0,db,w
       pi =4.*atan(1.)
       open(10,file='test')
       open(11,file='res')
       open(12,file='var')
c
c Eingabe der Parameter 
c
      read(10,*) na,nb,nc,dt
      read(10,*) m,b,v0,db
      w=2.*b*b

      do 1000 ia=1,na           
c
c alles nullsetzen
c
          do 100 i=1,nb
             do 100 j=1,3
                sig(j,i)=0.
                mu(j,i)=0.
100       continue
c     
c
c Anfangskonfigurationen
c
          print *,'Anfangskonfigurationen lesen',ia
          do 190 l2=1,nc
            print *,'Konfiguration lesen :',l2
            do 300 i=0,nb+1
               read(10,1010) (x(j,i,l2),j=1,3)
300         continue
c
c berechne Mittelwert und Varianz
c
            print *,'Mittelwerte + Var :'
            do 310 i=1,nb
               do 320 j=1,3
                  mu(j,i)=mu(j,i)+x(j,i,l2)/nc
                  sig(j,i)=sig(j,i)+x(j,i,l2)*x(j,i,l2)/nc
320            continue
310        continue
c
c Konfiguration abgearbeitet
c
190       continue

          print *,'Varianz korrigieren :',ia
c
c Varianz korrigieren
c
          sm=1.
          do 330 i=1,nb
             do 330 j=1,3
                sig(j,i)=sig(j,i)-mu(j,i)**2
                print *,sig(j,i)
                sm=sm*sqrt(2*pi*sig(j,i))
330       continue

c
c Testausgabe Mittelwert
          do 700 i=1,nb
             write(11,*) i,(mu(j,i),j=1,3)
             write(12,*) i,(sig(j,i),j=1,3)
700       continue

          goto 9999
c
c Konfiguration mit Massfunktional messen
c
          print *,'Messen Anfangskonf :',ia
          tm=0.
          do 200 l1=1,nc
             hm=0.
             print *,'Konfiguration :',l1
             do 210 i=1,nb 
               r=0.
               do 220 j=1,3
                  r=r+x(j,i,l1)**2
                  hm=hm
c freie Wirkung
     A            +m/db*(x(j,i+1,l1)-x(j,i,l1))**2          
c Massfunctional
c     B            -1/(2.*sig(j,i))*(x(j,i,l1)-mu(j,i))**2             
                   print *,i,j,hm
c                  print *,m/db*(x(j,i+1,l1)-x(j,i,l1))**2
c                  print *,-1/(2.*sig(j,i))*(x(j,i,l1)-mu(j,i))**2
220            continue
               r=sqrt(r)
c Potential 
               hm=hm+db*v(r)
210         continue
c 0 ter schritt
            r=0.
            do 230 j=1,3
               r=r+x(j,0,l1)**2
               hm=hm
     A           +m/db*(x(j,1,l1)-x(j,0,l1))**2                                       
230         continue
            r=sqrt(r)
            hm=hm+db*v(r)
            print *,'Wirkung minus Messfunktional :',hm
            tm=tm+exp(hm)
200       continue
c
c Messprozess fertig
c
          tm=tm/nc
          ids=ids+sm/tm
1000   continue
       ids=ids/na
c
c ids sollte jetzt das Pfadintegral gemessen ueber die
c Potentiale enthalten
c
c
c Normierungen:
c Pfadintegral und Potential
       ids=ids*(m/(2*pi*db))**(3./2.*nb)*(2*pi*b)**(3./2.)*v0**2
c Streulaenge
       ids=ids*((2*pi)**3/(4.*pi)*(2.*m))**2
c beta
       ids=ids*(2*pi*m/(nb*db))**(3./2.)
c
c Ausgabe 
c
       print *,'Streulaenge :'
       print *,sqrt(ids)
                       
        
1010   format(3e12.4)
9999   continue
       end
      

      real function v(r)
c
c potential (radialsymetrisch)
c
      real pi,m,b,v0,db,w 
      real r 
      common pi
      common/par/m,b,v0,db,w

      v=v0*exp(-r*r/w)
      end

      
