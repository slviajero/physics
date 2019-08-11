       program makep
       implicit none
c parameter:
       integer npolm,nlemax,nlamax
       real over
       parameter(npolm=10,nlemax=6,nlamax=7,over=9.)
C-----------------------------------------------------------
C Autor:   Stefan Lenz
C Datum:
C Version: 0.0
C-----------------------------------------------------------
C Dieses Program macht die Stuetzstellen fuer die
C numerische Integration der Vielfachstreugleichungen
C fuer zwei harmonische Oszillatoren und ein Target
C Es wird Gauss-Legendre und Gauss-Laguerre
C Quadratur zwischen den Polen bzw fuer grosse Impulse
C verwendet. Das Ausgabefile enthaelt die Stuetzstellen
C und Gewichte. An den physikalischen Polen wird
C Gewicht 0 ausgegeben.
C----------------------------------------------------------
C Parameter:
C npolm: maximale Zahl der physikalischen Pole
C nlemax: maximale Zahl der Legendre Punkte
C nlamax: maximale Zahl der Laguerre Punkte
C over: Differenz zwischen letztem Pol und hoechsten
c       p-Wert in Einheiten von pin
C
C lokales:
C npol: zahl der physikalischen pole
C pol: Polwerte
C nle: Zahl der Legendrepkt. zwischen zwei Polen
C nla: Zahl der Laguerrepkt.
c xle,wle: Stuetzstellen und Gewichte der Legendreploynome
c xla,wla: Stuetzstellen und Gewichte der Laguerrepolynome
C lambda : Reskalierungsfaktor 
       integer npol,nle,nla
       real pol(0:npolm)
       real xle(nlemax,nlemax)
       real wle(nlemax,nlemax)
       real xla(nlamax,nlamax)       
       real wla(nlamax,nlamax)
       real lambda
c
c triviales:
       real pi
       integer i,j
       real pmin,pmax,px,w

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

       common /param/mp,mt,hb,omega,pin,tau,nt,np

       pi = 4.*atan(1.)
c
c Wertetabelle:
c nach Abramowitz S916 und S923
c Legendre-Polynome
c  1.Polynom
        xle(1,1)=0.
        wle(1,1)=2.
c  2.Polynom
        xle(2,1)=-0.5774
        wle(2,1)= 1.0
        xle(2,2)= 0.5774
        wle(2,2)= 1.0
c  3.Polynom
        xle(3,1)=-0.7746
        wle(3,1)= 0.5556
        xle(3,2)= 0.0
        wle(3,2)= 0.8889
        xle(3,3)= 0.7746
        wle(3,3)= 0.5556
c  4.Polynom
        xle(4,1)=-0.8611
        wle(4,1)= 0.3479
        xle(4,2)=-0.3400
        wle(4,2)= 0.6521
        xle(4,3)= 0.3400
        wle(4,3)= 0.6521
        xle(4,4)= 0.8611
        wle(4,4)= 0.3479
c  5.Polynom
        xle(5,1)=-0.9062
        wle(5,1)= 0.2369
        xle(5,2)=-0.5385
        wle(5,2)= 0.4786
        xle(5,3)= 0.0000
        wle(5,3)= 0.5688
        xle(5,4)= 0.5385
        wle(5,4)= 0.4786
        xle(5,5)= 0.9062
        wle(5,5)= 0.2369
c   6.Polynom
        xle(6,1)=-0.9325
        wle(6,1)= 0.1713
        xle(6,2)=-0.6612
        wle(6,2)= 0.3608
        xle(6,3)=-0.2386
        wle(6,3)= 0.4679
        xle(6,4)= 0.2386
        wle(6,4)= 0.4679
        xle(6,5)= 0.6612
        wle(6,5)= 0.3608
        xle(6,6)= 0.9325
        wle(6,6)= 0.1713
c Laguerre-Polynome
c 2.Polynom:
        xla(2,1)= 0.5858
        wla(2,1)= 1.5333
        xla(2,2)= 3.4142
        wla(2,2)= 4.4510
c 3.Polynom
        xla(3,1)= 0.4158
        wla(3,1)= 1.0777
        xla(3,2)= 2.2943
        wla(3,2)= 2.7621
        xla(3,3)= 6.2900
        wla(3,3)= 5.6011
c 4.Polynom
        xla(4,1)= 0.3225
        wla(4,1)= 0.8327
        xla(4,2)= 1.7458
        wla(4,2)= 2.0481
        xla(4,3)= 4.5366
        wla(4,3)= 3.6311
        xla(4,4)= 9.3951
        wla(4,4)= 6.4871
c 5.Polynom
        xla(5,1)= 0.2636
        wla(5,1)= 0.6791
        xla(5,2)= 1.4134
        wla(5,2)= 1.6385
        xla(5,3)= 3.5964 
        wla(5,3)= 2.7694
        xla(5,4)= 7.0858
        wla(5,4)= 4.3157
        xla(5,5)=12.6408
        wla(5,5)= 7.2192
c 6.Polynom
        xla(6,1)= 0.2228
        wla(6,1)= 0.5735
        xla(6,2)= 1.1889
        wla(6,2)= 1.3693
        xla(6,3)= 2.9927
        wla(6,3)= 2.2607
        xla(6,4)= 5.7751
        wla(6,4)= 3.3505
        xla(6,5)= 9.8375
        wla(6,5)= 4.8868
        xla(6,6)=15.9829
        wla(6,6)= 7.8490
c 7. Polynom
        xla(7,1)= 0.1930
        wla(7,1)= 0.4965
        xla(7,2)= 1.0267
        wla(7,2)= 1.1776
        xla(7,3)= 2.5679
        wla(7,3)= 1.9182
        xla(7,4)= 4.9004
        wla(7,4)= 2.7718
        xla(7,5)= 8.1822
        wla(7,5)= 3.8412
        xla(7,6)=12.7342
        wla(7,6)= 5.3807
        xla(7,7)=19.3957
        wla(7,7)= 8.4054
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
c Ausgabefile oeffnen
c
       open(11,file='pdisc')
c
c Bestimme Zahle der Pole:
c
       npol = int(pin*pin/(2.*mp*hb*omega))+1
       print *,'Zahl der Pole',npol
       if (npol.gt.npolm) then
           print *,'Error: zu viele physikalische Pole:'
           print *,'npol  = ',npol
           print *,'npolm = ',npolm
           stop
       endif
c
c Berechne Impulswerte
c
       do 100 i=0,(npol-1)
          pol(i)=sqrt(pin*pin-2.*mp*hb*omega*i)
          print *,'Pol Nummer ',(i+1),' bei ',pol(i)
100    continue
       pol(npol)=0.
c
c Berechne Ordnung der Legendre und Laguerre
c Approximation
c
       nle=int((np-npol)/(npol+2.)*1.5)
       print *,'Zahl der Legendrepunkte :',nle
       if (nle.lt.1) then
          print *,'Error: Zahl der Legendrepkte kleiner Null'
          print *,'Zahl der Pole = ',npol
          stop
       endif
       if (nle.gt.nlemax) then
          print *,'Error: Zahl der Legendrepkte zu gross'
          print *,'Zahl der Pole = ',npol
          print *,'Zahl der Legendrepkte = ',nle
          stop
       endif
       nla=np-npol*(nle+1)
       print *,'Zahl der Laguerrepunkte ',nla
       if (nla.lt.1) then
          print *,'Error: Zahl der Laguerrepkte kleiner Null'
          print *,'Zahl der Pole           = ',npol
          print *,'Zahl der Legendrepkt    = ',nle
          stop
       endif
       if (nla.gt.nlamax) then
          print *,'Error: Zahl der Laguerrepkte zu gross'
          print *,'Zahl der Pole = ',npol
          print *,'Zahl der Legendrepkte = ',nle
          stop
       endif
c
c Berechne Integrationspunkte und -gewichte
c im Legendrebereich
c
       do 110 i=npol,1,-1
          pmin=pol(i)
          pmax=pol(i-1)
          do 120 j=1,nle
             px=(pmax-pmin)/2.*(xle(nle,j)+1)+pmin
             w=(pmax-pmin)/2.*wle(nle,j)
             write(11,*) px,w
120       continue             
          write(11,*) pmax,0.
110    continue            
c
c Berechne Intergrationspunkte und Gewichte im Laguerrebereich
c
c finde Reskalierungsfaktor
c (der hoechste Punkt liegt bei over*pin
       pmin=pol(0)
       lambda=over*pin/xla(nla,nla)
c
       do 130 i=1,nla
         px=pmin+xla(nla,i)*lambda
         w=wla(nla,i)*lambda
         write(11,*) px,w
130    continue
       close(11)
         
       end


