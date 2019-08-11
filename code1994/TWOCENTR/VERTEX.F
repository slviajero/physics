       program vertex
       implicit none
c parameter:
       integer ntm,npm
       parameter(ntm=10,npm=20)
c external functions:
       real horner
       external horner
C-----------------------------------------------------------
C Autor:   Stefan Lenz
C Datum:
C Version: 0.0
C-----------------------------------------------------------
C Dieses Programm berechnet die Vertexfunktionen fuer
C die Steuung an einem harmonischen Oszillator, der
C an das Projektil mit einer deltaartigen Tau-Matrix
C koppelt. Das hier sollte alle Infomationen ueber das
C Targetsystem enthalten.
C----------------------------------------------------------
C Parameter:
C ntm: maximale Zahl der Targetzustaende
C npm: maximale Zahl der Impulszustaende
c
c lokales:
c pwerte: enthaelt die Impulswerte die von makep
c         gemacht wurden
C pminus: enthaelt die Impulsdifferenzen
C pplus : enthaelt die Impulssummen
C eminus: Gaeusse an den entsprechenden pminus Werten
C eplus : Gaeusse an den entsprechenden pplus  Werten
C a     : Impulsskala des HO
C fak   : Fakultaeten bis nt
C c     : dimensionlose Normierungskonstanten des HO
C         sqrt(a/sqrt(pi)) ist schon rausgenommen
C d     : Polynomkoeffizienten
C fnp,fnm : Funktionswerte fuer plus und minus Argumente
C nk      : Zahl der Impulsdifferenzwerte
C
       real pwerte(npm)
       real pplus(npm*(npm+1)/2)
       real pminus(npm*(npm+1)/2)
       real eplus(npm*(npm+1)/2)
       real eminus(npm*(npm+1)/2)
       real a
       real fak(0:ntm)
       real c(0:ntm)
       real d(0:ntm)
       real fnp,fnm
       integer nk
C
C triviales:
       real void
       integer i,j,k
       integer n1,n2
       real pi
       real xm,xp

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

       pi=4.*atan(1.)
C
C liest die Parameter vom Eingabefile
C
       open(10,file='parin',ERR=1002)
       read(10,*,ERR=1003) mt,mp
       read(10,*,ERR=1004) hb
       read(10,*,ERR=1005) omega
       read(10,*,ERR=1006) pin
       read(10,*,ERR=1007) tau
       read(10,*,ERR=1008) nt,np
       close(10,ERR=1002)
       a=sqrt(mt*omega/hb)
       if (nt.gt.ntm) then
          print *,'Error: dieses Programm kann nur ',ntm
          print *,' Target-Basiszustaende'
          print *,' nt = ',nt
          stop
       endif
       if (np.gt.npm) then
          print *,'Error: dieses Programm kann nur ',npm
          print *,'Impuls-Basiszustaende'
          print *,' np = ',np
          stop
       endif
c
c Einlesen der Impulswerte
c
      open(11,file='pdisc',ERR=1000)
      do 100 i=1,np
         read(11,*,ERR=1000) pwerte(i),void
100   continue
      close(11,ERR=1000)
c
c Berechnung der Summen und Differenzen
c dividiert durch die entsprechenden Konstanten
C und der Exponetialfunktionen von den HO
C Wellenfunktionen
c
      k=1
      do 110 i=1,np
         do 120 j=1,i
            pplus(k)  = (pwerte(i)+pwerte(j))/(hb*a)
            pminus(k) = (pwerte(i)-pwerte(j))/(hb*a)
            eplus(k)  = exp(-(pplus(k)/2)**2)
            eminus(k) = exp(-(pminus(k)/2)**2)
            k=k+1
120      continue
110   continue
      nk=k-1
c
c Berechnung der Fakultaeten bis nt
c
      fak(0)=1.
      do 130 i=1,nt
         fak(i)=fak(i-1)*i
130   continue
      print *,fak
c
c Berechnung der HO Normierungskonstante
c
      do 140 i=0,(nt-1)
         c(i)=1/(sqrt(2**i*fak(i)))
140   continue
      print *,c
c
c Berechnung der Vertexfunktionen
c
      open(12,file='fvert',ERR=1001)
c Schleife mit n2>=n1      
      do 150 n1=0,(nt-1)
        n2=n1 
c
c Koeffizienten fuer n2=n1 machen
c
        do 170 j=0,n1
           d(j)=2.**(-j)*(-1)**j/(fak(n1-j)*fak(j)**2)
170     continue
160     continue
          print *,n1,n2,d
c
c Polynom fuer alle Impulswerte berechnen und speichern
c
          do 180 k=1,nk
             xp=pplus(k)
             xm=pminus(k)
             fnm=horner(xm,d,n1)*(-1)**n2*2**n1*xm**(n2-n1)
             fnm=fnm*eminus(k)*c(n1)*c(n2)
             fnp=horner(xp,d,n1)*(-1)**n2*2**n1*xp**(n2-n1)
             fnp=fnp*eplus(k)*c(n1)*c(n2)
             write(12,*,ERR=1001) fnm,fnp
180       continue
c
c
c Koeffizienten fuer n2+1 machen
c
        if (n2.ge.(nt-1)) goto 200
           n2=n2+1
           do 190 j=0,n1
               d(j)=d(j)/(n2-n1+1+j)
190        continue
c
c Ende der n2-Schleife
c
        goto 160
200     continue
c
c Ende der n1 Schleife
c
150   continue
      close(12,ERR=1001)
      goto 9999
c
c io error traps
c
1000  continue 
         print *,'Error : Lesefehler bei Einlesen der Impulse'
         print *,'Impulswert No.: ',i
         stop
1001  continue 
         print *,'Error : Schreibfehler '
         stop
1002  continue   
         print *,'Error : parin File '
         stop
1003  continue   
         print *,'Error : parin File: Massen '
         stop
1004  continue   
         print *,'Error : parin File: hb '
         stop
1005  continue   
         print *,'Error : parin File: Omega '
         stop
1006  continue   
         print *,'Error : parin File: pin '
         stop
1007  continue   
         print *,'Error : parin File: tau '
         stop
1008  continue   
         print *,'Error : parin File: nt,np '
         stop
9999  continue
      end

      real function horner(x,a,n)
      integer n
      real x
      real a(0:n)
C
C rechnet Polynome nach dem Hornerschema aus
C a: Feld mit Polynomkoeffizienten
C n: Grad des Polynoms
C x: Argument
C lokales:
      integer i
      real y

      y=a(n)
      if (n.gt.0) then
        do 100 i=(n-1),0,-1
          y=y*x+a(i)
100     continue
      endif

      horner=y
      end



