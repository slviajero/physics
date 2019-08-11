       program fold
       implicit none
c parameter:
       integer nconf,nzent
       parameter(nconf=50000,nzent=2)
C
C nconf: Zahl der Konfigurationen
C nzent: Zahl der Zentren 
C-----------------------------------------------------------
C Autor:   Stefan Lenz
C Datum:
C Version: 0.0
C-----------------------------------------------------------
C Loest das Foldymodell in zwei Raumdimensionen
C (und gibt 2te Ordnung Born aus)
C-----------------------------------------------------------
C
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
C
C Definitionen des Zufallszahlengenerators
C
       DOUBLE PRECISION X,XK
       COMPLEX TT
c
c lokales:
C field : Konfiguration
C rij   : Abstandsmatrix
C mg    : Streumatrix
C incv  : einlaufende ebene Welle
C outvf : auslaufende Welle Vorwaertsstreuung
C outvb : auslaufende Welle Rueckwaertsstreuung
C nc    : Konfigurationszaehler
C b     : Loesungsvektor
C h     : Hilfsgroessen fuer die Loesungsroutine 
C tf,tb : T Matrix in Vorwaerts und Rueckwaertsrichtung 
C t,r   : Transmission+Reflexionsamplituden
C Entsprechen t1... , t2... fuer Bornterme
C
       real field(nzent)
       real rij(nzent,nzent)
       complex mg(nzent,nzent)
       complex incv(nzent)
       complex outvf(nzent)
       complex outvb(nzent)
       integer nc
       complex b(nzent),h(nzent),hs
       complex t,r,tb,tf
       complex tb1,tf1
       complex tb2,tf2
       character*(100) line
C
C Funktionen:
        real ran,rgau
        external ran,rgau

C
C triviales:
       integer i,j,k
       real pi
       complex ic

       COMMON/RG/X,XK
       common /param/mp,mt,hb,omega,pin,tau,nt,np,pi

       pi = 4.*atan(1.)
       ic = (0. , 1.)
C 
C Initialisierung des Zufallszahlengenerators
C
      X=PI*1.D11
      XK=0.D00
      DO 10 K=1,20
        TT=RAN()
10    CONTINUE
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
c Resultatfile umkopieren
       open(21,file='fold.res')
       open(22,file='fold.res.new')
 28    read(21,8000,END=29) line
       write(22,8000) line
       goto 28
 29    continue
c
c einiges Nullsetzen 
       tb=(0.,0.)
       tf=(0.,0.)
       tb1=(0.,0.)
       tb2=(0.,0.)
       tf1=(0.,0.)
       tf2=(0.,0.)
c
c Hauptschleife
       do 1000 nc=1,nconf
c
c ziehe Konfiguration
          do 100 i=1,nzent
             field(i)=rgau()*sqrt(hb/(2.*mt*omega))
100       continue
c
c berechne Abstandsmatrix
          do 110 i=1,nzent-1
             do 110 j=i+1,nzent
                rij(i,j)=abs(field(i)-field(j))
                rij(j,i)=rij(i,j)
110       continue
c
c berechne Greensfunktionsmatrix
          do 120 i=1,nzent-1
             do 120 j=i+1,nzent
                mg(i,j)=tau*ic*mp/(hb*pin)*exp(ic*pin*rij(i,j))
                mg(j,i)=mg(i,j)
120       continue
          do 130 i=1,nzent
             mg(i,i)=(1.,0.)
130       continue
c
c berechne einlaufenden Vektor
          do 140 i=1,nzent
             incv(i)=exp(ic*pin*field(i))
140       continue
c
c berechne die auslaufenden Vektoren
          do 150 i=1,nzent
             outvf(i)=exp(-ic*pin*field(i))
             outvb(i)=exp(ic*pin*field(i))
150       continue
c
c hier werden die erste und zweite Ordnung Bornapproximation
c berechnet
c erste Ordnung:
	  do 160 i=1,nzent
	     tb1=outvb(i)*incv(i)+tb1
	     tf1=outvf(i)*incv(i)+tf1
160       continue
c
c zweite Ordnung:
          do 170 i=1,nzent
	     h(i)=(0.,0.)
	     do 180 j=1,nzent
		if (i.ne.j) then
		   h(i)=-(mg(i,j)*incv(j))+h(i)
                endif
 180         continue
             tb2=tb2+h(i)*outvb(i)
	     tf2=tf2+h(i)*outvf(i)
170       continue
c
c loese Gleichungssystem: 
C  Methode : Bronstein, Taschenbuch der Mathematik
C           S.735 ( Gaussches Eliminationsverfahren )
c
C Transform A to triangular shape
          do 200 i=1,(nzent-1)
             do 230 j=(i+1),nzent
                h(j)=-mg(j,i)/mg(i,i)
                do 220 k=1,nzent
                  if (k.eq.i) then
                     mg(j,k)=0.0
                  else
                     mg(j,k)=mg(j,k)+h(j)*mg(i,k)
                  endif
 220            continue
                incv(j)=incv(j)+h(j)*incv(i)
 230         continue
 200      continue
C Generate solution vector
          b(nzent)=incv(nzent)/mg(nzent,nzent)
          do 300 i=(nzent-1),1,-1
             hs=0.0
             do 310 j=1,(nzent-i)
                hs=hs+mg(i,i+j)*b(i+j)
 310         continue
          b(i)=(incv(i)-hs)/mg(i,i)
 300      continue
c
c Loesungsvektor auf auslaufende Welle projizieren
          hs=(0.,0.)
          do 400 i=1,nzent
             hs=hs+b(i)*outvf(i)
400       continue 
          tf=tf+hs
          hs=(0.,0.)
          do 410 i=1,nzent
             hs=hs+b(i)*outvb(i)
410       continue
          tb=tb+hs
1000   continue
       tb=tau*tb/nconf
       tf=tau*tf/nconf 
       tb1=tau*tb1/nconf
       tf1=tau*tf1/nconf 
       tb2=tau*tb2/nconf
       tf2=tau*tf2/nconf 
c
c jetzt stehen in tb und tf die mittleren T matrizen
c
       r = -mp*ic/(hb*pin)*tb
       t = (1.,0.)-(mp*ic)/(hb*pin)*tf
       write(22,8007) tau,pin,t,r

       print *,'Exaktes Resultat  :'
       print *,'Reflexionsamplitude: ',r,abs(r)**2
       print *,'Transmissionsamplitude: ',t,abs(t)**2
       print *,'T+R = ',abs(r)**2+abs(t)**2

       r = -mp*ic/(hb*pin)*tb1
       t = (1.,0.)-(mp*ic)/(hb*pin)*tf1
       print *,'Erste Ordnung Born:'
       print *,'Reflexionsamplitude: ',r,abs(r)**2
       print *,'Transmissionsamplitude: ',t,abs(t)**2
       print *,'T+R = ',abs(r)**2+abs(t)**2
 
       r = -mp*ic/(hb*pin)*(tb2)
       t = -(mp*ic)/(hb*pin)*(tf2)
       print *,(tb2+tf2)/2,(tf2-tb2)/2
       print *,'Zweite Ordnung Born:'
       print *,'Reflexionsamplitude: ',r,abs(r)**2
       print *,'Transmissionsamplitude: ',t,abs(t)*2
       print *,'T+R = ',abs(r)**2+abs(t)**2
c
c Ergebnissfiles umkopieren
c
        close(21)
        close(22)
        open(21,file='fold.res.new')
        open(22,file='fold.res')
 2000   read(21,8000,END=2001) line
        write(22,8000) line
        goto 2000
 2001   continue

c
c Formats
c
 8000  format(a100)
 8007  format(7e12.4)
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
        rgau=sqrt(-2.*alog(ran()))*cos(2*3.142*ran())
      end



