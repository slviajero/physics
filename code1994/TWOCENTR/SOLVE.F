       PROGRAM SOLVE
       implicit none
c parameter:
       integer ntm,npm
       parameter(ntm=1,npm=20 )
c external functions:
C-----------------------------------------------------------
C Autor:   Stefan Lenz
C Datum:
C Version: 0.0
C-----------------------------------------------------------
C Dieses Programm loest die Vielfachstreugleichungen
C in einem eingeschraekten Hilbertraum durch  Trigonalisieren
C der Rueckstreumatrix. Es benoetigt als Eingabe die Files
C pdisc und fvert, die von den Programmen makep und vertex 
C erstellt werden. Wegen des hohen Speicherplatzbedarfs
C sollte immer genau die Anzahl Basisvektoren verwendet
C werden,  die dimensioniert wurde.
C----------------------------------------------------------
C Speicherplatzbedarf: es werden (ntm**4)*(npm*2)
C Speicherworte benoetigt und die Matrix aufzunehmen
C 
C----------------------------------------------------------
C Ein Wort zu den Indexdateien ki und nvi:
C Nachdem in der Vertexdatei die Vertizes in
C gepackter Form stehen, (dh.: die Symmetrie
C unter Vertauschung von n und np, sowie 
C die Symmetrieeigenschaften unter Vertauschung
C con k,kp werden ausgenutzt) ist die Selektion
C der richtigen Vertizes tricky. Daher wurden 
C diese Indexdateien eingefuehrt. Aenderungen 
C hier sollten nur mit grosser Vorsicht 
C angebracht werden
C
C Parameter:
C ntm: maximale Zahl der Targetzustaende
C npm: maximale Zahl der Impulszustaende
c
c lokales:
c pwerte: enthaelt die Impulswerte die von makep
c         gemacht wurden
C nk    : Zahl der Impulsdifferenzwerte
C fnp   : Vertexfkt fuer P+ 
C fnm   : Vertexfkt fuer P-
C npol  : HO Quantenzahl bei der Pole existieren
C pole  : Enthaelt die Indizes der Pol-pwerte als Fkt von n1+n2
C rsm   : Rueckstreumatrix
C incv  : einlaufende Welle
C outv  : auslaufende Welle
C h     : Hilfsvektor fuer Trigonalisierung
C sol   : Loesungsvektor
C prin,prout : einlaufende und auslaufende Targetparitaet
C pr1    : Pariaet von Target 1
C pr2    : Paritaet von Target 2
C pr12p  : Paritaet von Zwischenzustand
C gp     : Gesammtparitaet
C al     : einlaudende Projetilparitaet 
C alp    : einlaufende Projetilparitaet
C alpp   : Zwischenzustand-Projetilparitaet
C s1p - s2m : Selektionsfaktoren fuer die
C             Uebergaenge (S47.)
C g      : Greensfunktion in der diskretisierten
C          Approximation
C nvi    : Indexdatei zur Auffindung der Vertizes (n1,n1p)
C ki     : Indexdatei zur Auffindung der Vertizes (k,kp)
C irsm   : Indexdatei zum Auffinden eines Elements aus rsm
C ikin   : Index auf elastischen Pol
C fpar   : enthaelt (-1)**n
C fipar  : enthaelt ic**(n)
C tplus,tminus : Loesungen fuer Paritaet plus und minus
C tplus1,tminus1 : Born Approx. 1te Ordnung
C tplus2,tminus2 : Born Approx. 2te Ordnung
C t,r          : Refexions und Transmissionskoeffizient
C tpm,tmm : Matrixelemente fuer die Testroutine
C pfaktor : enthaelt die richtigen Faktoren fuer
C           die Paritaetszerlegungen der t Matrizen
c line    : Eingabestring
       real pwerte(npm)
       real w(npm)
       integer nk
       real fnp(npm*(npm+1)/2,0:ntm*(ntm+1)/2-1)
       real fnm(npm*(npm+1)/2,0:ntm*(ntm+1)/2-1)
       integer npol
       integer pole(0:npm-1)
       complex rsm(ntm*ntm*npm,ntm*ntm*npm)
       complex incv(ntm*ntm*npm)
       complex outv(ntm*ntm*npm)
       complex h(ntm*ntm*npm)
       complex sol(ntm*ntm*npm)
       integer pr1,pr2,pr12p,prin,prout
       integer gp
       integer al,alp,alpp
       integer s1p,s1m,s2p,s2m  
       complex g(0:(ntm-1),0:(ntm-1),npm)
       integer nvi(0:(ntm-1),0:(ntm-1))
       integer ki(npm,npm)
       integer irsm(npm,0:(ntm-1),0:(ntm-1))
       integer ikin
       integer fpar(0:ntm*ntm+1)
       complex fipar(0:ntm*ntm+1)
       complex tplus,tminus
       complex tplus1,tminus1
       complex tplus2,tminus2
       complex t,r
       complex tpm,tmm
       complex pfaktor(-1:1,-1:1)
       character*(100) line
C
C triviales:
       integer i,j,k
       integer n1,n2,nn
       integer ni,nip
       integer n1p,n2p
       integer k1
       integer k1p
       integer kpp
       integer n1v,n2v
       integer ppv
       integer psig
       real pi
       complex gg
       complex gg1,t1,gg2,t2
       complex ic
       complex ir1,ir2 
       integer pnuml
       complex hc

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

c
c Allgemeine Konstanten:
C
       pi=4.*atan(1.)
       ic=(0.,1.)
C      
C Paritaeten
       fpar(0)  =1
       fipar(0) =(1.,0.)
       do 10 i=1,ntm*ntm+1
	   fpar(i)  =  fpar(i-1)*(-1)
	   fipar(i) =  fipar(i-1)*ic
10     continue
c
c Paritaetsfaktoren vorbesetzen
c
        pfaktor(-1,-1)= 1.
        pfaktor( 1,-1)=-ic
        pfaktor(-1, 1)= ic
        pfaktor( 1, 1)= 1.
C
C oeffne Hilfdatei
c
        open(20,file='ftn20')
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
c oeffne Ergebnissdatei
c und kopiere um
       open(21,file='sol.res')
       open(22,file='sol.res.new')
 28    read (21,8000,END=29) line
       write(22,8000) line
       goto 28
 29    continue
c
c Aufbau der rsm Indexdatei
c
      nn=1
      do 30 n1=0,nt-1
	 do 31 n2=0,nt-1
	    do 32 k=1,np
	       irsm(k,n1,n2)=nn
	       nn=nn+1
32          continue
31       continue
30    continue
c
c Einlesen der Impulswerte
c
      open(11,file='pdisc',ERR=1000)
      do 100 i=1,np
	 read(11,*,ERR=1000) pwerte(i),w(i)
100   continue
      close(11,ERR=1000)
c      
c  Anzahl der Impulswerte
c
      nk=np*(np+1)/2
c
c Ermitteln der Polwerte
c hier werden die Energien bei denen
c Pole existieren koennen festgestellt
c
      npol=int(pin**2/(2*mp*hb*omega))
      n1=0
      do 89 i=np,1,-1
	  if (w(i).eq.0.) then 
	     pole(n1)=i
	     n1=n1+1
	  endif
89    continue
      if ((n1-1).ne.npol) then
	 print *,'Error : Pole falsch gezaehlt'
	 stop
      endif
c
c Aufstellen der Greensfunktion
c-----------------------------------------------
c der Faktor 2 bei w(k) kommt von der Trafo
c Integral -unendlich bis +unendlich auf
C Integral 0 bis unendlich
c-----------------------------------------------
c 1. Besetzen der nichtsingulaeren Stellen
      do 90 n1=0,(nt-1)
	 do 91 n2=0,(nt-1)
	   do 92 k=1,np
	      if (w(k).ne.0.) then
		 g(n1,n2,k)=cmplx(
     A         2.*w(k)/((pin**2-pwerte(k)**2)/(2.*mp)-(n1+n2)*hb*omega)
     B                           )
	      else
		 g(n1,n2,k)=(0.,0.)
	      endif
92         continue
91       continue
90    continue
C 2. Korrekturterme einfuegen
      do 93 i=0,npol
	 do 94 n2=0,i
	    n1=i-n2
	    gg=(0. , 0.)
	    do 95 j=1,np
	       gg=gg+g(n1,n2,j)
95          continue
	    gg=-gg-2.*ic*pi*mp/(pwerte(pole(i)))
	    g(n1,n2,pole(i))=gg
94       continue
93    continue
c
c 
c      print *,'Impulswerte:'
c      print *,pwerte
c      print *,'Greensfunktionsmatrix:'
c      print *,g
      ikin = pole(0)
c
c Einlesen der Vertexfunktionen
c hier wird auch gleichzeitig die Indexdatei 
C nvi aufgebaut. Man muss dann nur noch 
C den entsprechenden Index aus nvi abholen
c      
      nn=0
      open(12,file='fvert',ERR=1001)
      do 150 n1=0,(nt-1)
	 do 160 n2=n1,(nt-1)
	   do 170 k=1,nk
	     read(12,*,ERR=1001) fnm(k,nn),fnp(k,nn)
170        continue
	   nvi(n1,n2)=nn
	   nvi(n2,n1)=nn
	   nn=nn+1
160      continue
150   continue
      close(12,ERR=1001)
c
c Aufbau der Indexdatei fuer k,kp
c
	nn=1  
	do 151 k=1,np
	   do 152 j=1,k
	      ki(k,j)=nn
	      ki(j,k)=nn
	      nn=nn+1
152        continue
151     continue
c
c Paritaetsschleife
c
       do 20 gp=-1,1,2

c
c Aufstellen der Rueckstreumatrix
c
	do 500 n1=0,(nt-1)
	do 500 n2=0,(nt-1)
	do 500 n1p=0,(nt-1)
	do 500 n2p=0,(nt-1)

c
c erzeuge Index fuer Targetfreiheitsgrade
c das sagt welcher k,kp Block gerade 
C bearbeitet wird
c
	ni  = irsm(1,n1,n2)-1
	nip = irsm(1,n1p,n2p)-1
c
c erzeuge Index fuer Vertexfunktion
c das sagt aus welchem Teil des Vertexarray
c die Funktionen genommen werden sollen
c (der groessere Index wird zuert hochgezaehlt!)
c
       
	n1v=nvi(n1p,n1)
	n2v=nvi(n2p,n2)
    
c
c ermittele Paritaeten fuer die einzelnen Kombinationen
c von Target und Projetilfkt.
c sowie die Selektionsfaktoren, die entscheiden
c ob fnm oder fnp fuer den Uebergang verwendet 
c wird
c
	prin  = fpar(n1+n2)
	prout = fpar(n1p+n2p)
	pr1   = fpar(n1+n1p)
	ir1   = fipar(n1+n1p)        
	pr12p = fpar(n1+n2p)
	pr2   = fpar(n2+n2p)
	ir2   = fipar(n2+n2p)
	al     = gp*pr1
	alp    = gp*pr2
	alpp   = gp*pr12p
	s1p    = (1+pr1*alp*alpp)
	s1m    = (alp+pr1*alpp)
	s2p    = (1+pr2*al*alpp)
	s2m    = (al+pr2*alpp)
c
c Aufstellen der k-Matrix 
c
	do 300 k1=1,np
	   do 300 k1p=1,np
c
c linke Greensfkt:
	   gg1=g(n1p,n2p,k1p)
c
c intermediaere Integration
	   gg=(0.,0.)
	   do 310 kpp=1,np 
c
c mache Indizes fuer das auslesen der k-werte aus
c fnm,fnp
c das ist so geordnet, dass der groessere Wert
c langsamer und der groessere schneller laueft
c psig ist das Vorzeichen, dass von den negative
c Impulsdifferenzen herkommt, es wird aus den
C Impulswerten und der Paritaet von n1,n1p bestimmt.
C 
C hier ist nochwas mit der Ausgangsparitaet verkehrt!!
	      if (kpp.ge.k1p) then 
		 psig=+1
	      else
		 psig=pr1
	      endif
	      ppv=ki(k1p,kpp)
	      t1=(fnp(ppv,n1v)*s1p+fnm(ppv,n1v)*psig*s1m)
     A           *pfaktor(alpp,alp)
	      gg2=g(n1,n2p,kpp)
	      if (kpp.gt.k1) then 
		 psig=pr2
	      else
		 psig=+1
	      endif
	      ppv=ki(k1,kpp)
	      t2=(fnp(ppv,n2v)*s2p+fnm(ppv,n2v)*psig*s2m)
     A           *pfaktor(al,alpp)
	      gg=gg+t1*gg2*t2
c
c Test der onshell Matrixelemente
c
c              if ((k1.eq.ikin).and.(k1p.eq.ikin)) then
c                 print *,'Impuls : ',kpp,pwerte(kpp)
c                 print *,'Mat:',t1*gg2*t2*(tau/(8.*pi))**2
c                 print *,'Paritaeten :'
c                 print *,al,alp,alpp
c		 print *,'Selektionsfaktoren :'
c		 print *,' 1: ',s1p,s1m
c		 print *,' 2: ',s2p,s2m
c		 print *,'Ts :'
c		 print *,t1*tau/(8*pi),t2*gg2*tau/(8*pi)
c		 print *,pr1,pr2
c              endif 
310        continue
c
c restliche Faktoren drauf  (siehe S47)
c
           gg=gg*alpp*gg1*(tau/(8.*pi))**2/fipar((2-al-alp)/2)
c
c Testroutine: 
c einzelne Kanael koennen hier rausgenommen werden.
c dazu ist zu sagen: die Differenz von 1ter und zweiter
c Ordung Bornapproximation ist genau das entsprechene
c Matrixelment in der Rueckstreumatrix, falls man das entprechende
c Bein amputiert, daher
            if (((n1+n2+n1p+n2p).eq.0).and.(k1p.eq.ikin)
     A             .and.(k1.eq.ikin)) then
               if (gp.eq.1) then 
                  tpm=2*pi*gg/gg1
               else
                  tmm=2*pi*gg/gg1
               endif
            endif 
c
c Eintragen des Matrixelements
c
	   if ((nip+k1p).eq.(ni+k1)) then
	         rsm(nip+k1p,ni+k1)=(1.,0.)-gg
           else
		 rsm(nip+k1p,ni+k1)=-gg
           endif
300     continue
500     continue
c       write(20,*) rsm
c
c an dieser Stelle ist die Matrix erstellt,
c jetzt kommen die Vektoren an die Reihe
c
c
c zuerst wird der einlaufende vektor aufgestellt, das
c Projetktil trifft zunaechst Zentrum 2 oder porpagiert von
c Zentrum 1 zu 2, es befindet sich im kin Kanal
c elastischer Kanal:
c einlaufende Welle bei 2:
c besetze alle vor
	do 610 i=1,np*nt*nt
	    incv(i)=(0.,0.)
610     continue
c pole(0) = niedrigsten Polwert = elastischer Kanal
	incv(ikin)=(1.0,0.)     
c
c einlaufende Welle bei 1
	n2=0
C
C Projetil Eingangspartitaet
	al=gp
	do 620 n1=0,nt-1
c
c Projetil Ausgansparitaet
	   alp=fpar(n1)*gp
	   do 630 k=1,np
C
C Index in Vektor setzen
	      ni=irsm(k,n1,n2)
C
C Vertex selektieren
	      ppv=ki(k,ikin)
	      n1v=nvi(n1,0)
C
C Vorzeichen fuer Vertex generieren
	      if (k.lt.ikin) then
		  psig=fpar(n1)
	      else
		  psig=1
	      endif
c
c Vektorelement besetzen
	      incv(ni)=incv(ni)+g(n1,n2,k)
     A                 *(tau/(8.*pi))/fipar((2-alp-al)/2)
     B                 *2.*(fnp(ppv,n1v)+fnm(ppv,n1v)*alp*psig)
     C                 *fipar(n1)*pfaktor(al,alp)
630        continue
620     continue
c
c auslaufende Welle bei 2:
C (Etwas hastig gemacht, kann gut falsch sein)
c besetze alle vor
	do 710 i=1,np*nt*nt
	    outv(i)=(0.,0.)
710     continue
	n1=0
C
C Projetil Ausgangspartitaet
	al=gp
	do 720 n2=0,nt-1
c
c Projetil Eingangsparitaet
	   alp=fpar(n2)*gp
	   do 730 k=1,np
C
C Index in Vektor setzen
	      ni=irsm(k,n1,n2)
C
C Vertex selektieren
	      ppv=ki(k,ikin)
	      n1v=nvi(n2,0)
C
C Vorzeichen fuer Vertex generieren
	      if (k.lt.ikin) then
		  psig=fpar(n2)
	      else
		  psig=1
	      endif
c
c Vektorelement besetzen
	      outv(ni)=outv(ni)+
     A                  (tau/(8.*pi))/fipar((2-alp-al)/2)
     B                 *2.*(fnp(ppv,n1v)+fnm(ppv,n1v)*gp*psig)
     C                 *fipar(n2)
730        continue
720     continue
c 	print *,'Paritaet :',gp
c	print *,outv
c       print *,'Matrix '
c       print *,rsm
c        print *,'Inhomog. Vektor '
c        print *,incv
c
c Born-Tests:
c Beachte: die ersten beiden Ordnungen Bornapproximation
c          testen nur die Vektoren und nicht die Matrix
c
c Berechnung der 1.ten Ordnung:
c Der zusaetzlich Faktor 2*pi kommt davon, dass ein
c und auslaufende Welle auf 2*pi und nicht auf eins
c normiert sind
        if (gp.eq.1) then
	   tplus1=2*pi*outv(ikin)  
        else
	   tminus1=2*pi*outv(ikin)
        endif
c
c Berechnung der Zweiten Ordnung
C (inklusive erster!!)
c
        hc=(0.,0.)
	do 750 i=1,nt*nt*np
            hc=hc+incv(i)*outv(i)
            if (i.lt.ikin) then
c                  print *,'Dir:',incv(i)*outv(i)
            elseif (i.eq.ikin) then
c                  print *,'Dir:',(incv(i)-1.)*outv(i)
            else
c                  print *,'Dir:',incv(i)*outv(i)
            endif
750     continue
        if (gp.eq.1) then
	   tplus2=2*pi*hc  
        else
	   tminus2=2*pi*hc
        endif
c
C Loesung des Gleichungssystems:
C dieser Programmteil wurde aus pn entnommen
C
C  Solve linear equations 
C

C
C  Auxilliary variables
      pnuml=nt*nt*np
C
C  Solve linear equations
C  Solution : b
C  Method : Bronstein, Taschenbuch der Mathematik
C           S.735 ( Gaussches Eliminationsverfahren )
c
C Transform A to triangular shape
	  do 1200 i=1,(pnuml-1)
	     do 1230 j=(i+1),pnuml
		h(j)=-rsm(j,i)/rsm(i,i)
		do 1220 k=1,pnuml
		  if (k.eq.i) then
		     rsm(j,k)=0.0
		  else
		     rsm(j,k)=rsm(j,k)+h(j)*rsm(i,k)
		  endif
 1220            continue
		incv(j)=incv(j)+h(j)*incv(i)
 1230         continue
 1200      continue
C Generate solution vector
	  sol(pnuml)=incv(pnuml)/rsm(pnuml,pnuml)
	  do 1300 i=(pnuml-1),1,-1
	     hc=(0.,0.)
	     do 1310 j=1,(pnuml-i)
		hc=hc+rsm(i,i+j)*sol(i+j)
 1310        continue
	  sol(i)=(incv(i)-hc)/rsm(i,i)
 1300     continue
c
c jetzt sollte in sol die exakte Loesung des Problems stehen.
C das Ergebnis muss jetzut noch auf eine auslaufende Welle
C projiziert werden
c
	   hc=(0.,0.)
	   do 1400 i=1,pnuml
	      hc=hc+sol(i)*outv(i)
1400       continue
c
c hc enthaelt jetzt das T-Matrixelement fuer Paritaet gp
c Der zusaetzlich Faktor 2*pi kommt davon, dass ein
c und auslaufende Welle auf 2*pi und nicht auf eins
c normiert sind
	   if (gp.eq.1) then
		tplus=2*pi*hc
	   else
		tminus=2*pi*hc
	   endif
c
c Ende der Schleife gp
c
20      continue
c
c Berechnung von Reflexions und Transmissionskoefizient
c
	t=(1.,0.)-mt*ic/(hb*pin)*2.*(tplus+tminus)
	r=-mt*ic/(hb*pin)*2.*(tplus-tminus)
        write(22,8007) tau,pin,t,r

	print *,'Exaktes Resultat:'
	print *,tplus,tminus
	print *,'Reflexionsamplitude: ',r,abs(r)**2
	print *,'Transmissionsamplitude :',t,abs(t)**2
	print *,abs(r)**2+abs(t)**2

	t=(1.,0.)-mt*ic/(hb*pin)*2.*(tplus1+tminus1)
	r=-mt*ic/(hb*pin)*2.*(tplus1-tminus1)
       
	print *,'Born, 1. Ordnung :'
	print *,tplus1,tminus1
	print *,'Reflexionsamplitude: ',r,abs(r)**2
	print *,'Transmissionsamplitude :',t,abs(t)**2
	print *,abs(r)**2+abs(t)**2

	t=-mt*ic/(hb*pin)*2.*(tplus2+tminus2-tplus1-tminus1)
	r=-mt*ic/(hb*pin)*2.*(tplus2-tminus2-tplus1+tminus1)
       
         print *,'Born, 2te Ordnung :' 
         print *,2*(tplus2-tplus1),2*(tminus2-tminus1)
         print *,'Reflexionsamplitude: ',r,abs(r)**2
         print *,'Transmissionsamplitude :',t,abs(t)**2
         print *,abs(r)**2+abs(t)**2

c	t=-mt*ic/(hb*pin)*2.*(tpm+tmm)
c	r=-mt*ic/(hb*pin)*2.*(tpm-tmm)
       
c        print *,'Born, 2te Ordnung aus Matrix:' 
c        print *,tpm,tmm
c        print *,'Reflexionsamplitude: ',r,abs(r)**2
c        print *,'Transmissionsamplitude :',t,abs(t)**2
c        print *,abs(r)**2+abs(t)**2
	goto 9999
c
c io error traps
c
1000  continue 
	 print *,'Error : Lesefehler bei Einlesen der Impulse'
	 print *,'Impulswert No.: ',i
	 stop
1001  continue 
	 print *,'Error : fvert Lesefehler '
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
c
c Ergebnissfiles umkopieren
c
      close(21)
      close(22)
      open(21,file='sol.res.new')
      open(22,file='sol.res')
 2000 read(21,8000,END=2001) line
      write(22,8000) line
      goto 2000
 2001 continue
c
c Formats
c
 8000 format(a90)
 8002 format(2e12.4)
 8003 format(3e12.4)
 8004 format(4e12.4)
 8005 format(5e12.4)
 8006 format(6e12.4)
 8007 format(7e12.4)

      end

   

