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
c weiter modifikation: streuung an a zentren
c-----------------------------------------------------------
       implicit none 
c arraydimensionen
       integer nbmax,amax
       parameter(nbmax=8000,amax=20)
c allegmeine parameter
       real pi
c externe Funktionen
       real ran,rgau,v,vf
       external ran,rgau,v,vf
c diverse summen
       real tt,st,sv,svv,q2av,xq2av,q2avi,xq2avi
c simulationsgroessen
       real db
       real ds,dt
       integer nc,n0,nb
       integer ng
c schleifen
       integer i,k1,j
       integer l1,l2,la
c Parameter des Projectils:
       real m,b,v0,w,mb
c Parameter des Oszi:
       real omega,mo,xo,mob,f1,f2
       integer a
       real alpha
c Anfangskonfiguration:
       real xi(3),xf(3)
       real qi(3,amax),qf(3,amax)
c diverse flags
       logical cwrite
       logical twrite
       logical hwrite
c
c Felder fuer Langevin 
c      
       real r2(0:nbmax+1,amax), r2t(0:nbmax+1) 
       real x(3,0:nbmax+1), dx(3,0:nbmax+1)
       real q(3,0:nbmax+1,amax), dq(3,0:nbmax+1,amax)       
       real dw(3,0:nbmax+1)
       real dwq(3,0:nbmax+1,amax)
       real z(3)
       real zq(3,amax)
c
c Ergbnisse
c
c direkt:
       real savd,sad
c indirekt:
       real savid,said
c endergebnis: Verhältnis des Wirkungsquerschnitts von vollem Problem zu efektivem Potential
       real sr
c
c Zufallsgenerator 
c
       double precision xr,xk
       logical buildin
c
c globale variablen 
c
       COMMON/RG/Xr,XK,BUILDIN
       common/const/pi
       common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2,a,alpha
c
c parameter
c
       pi =4.*atan(1.)
c
c testflags setzen
c
       cwrite=.false.
       twrite=.true.
       hwrite=.false.
c
c Random Number flag
c
	buildin=.true.

c
c files
c
       open(9,file='result')
       open(10,file='testd')
       open(11,file='testid')
       open(12,file='error')
       open(14,file='ratio')
       open(15,file='halfconfig')
       open(16,file='radius')
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
      if (nb.lt.1) then
	  print *,'nb<1 !'
	  goto 9999
      elseif (nb.gt.nbmax) then
	  print *,'nbmax=',nbmax
	  goto 9999
      endif
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
      print *,'Zahl der Zentren :'
      read *,a
      if (a.lt.1) then
	  print *,'a<1 !'
	  goto 9999
      elseif (a.gt.amax) then
	  print *,'amax=',amax
	  goto 9999
      endif


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
      
c Protokollfile schreiben
      write(9,*) "<SIMULATION>"
      write(9,*) "Lang2hoa: Version $Revision: 1.17 $ "
      write(9,*) "<PAR>"
      write(9,*) "nbeta ", nb
      write(9,*) "dbeta ", db
      write(9,*) "b ", b
      write(9,*) "v0 ", v0
      write(9,*) "m ", m 
      write(9,*) "mo ", mo 
      write(9,*) "omega ", omega 
      write(9,*) "a ",a
      write(9,*) "nc ",nc
      write(9,*) "n0 ",n0
      write(9,*) "dt ",dt
      write(9,*) "</PAR>"
      write(9,*) "<RESULT>"
      write(9,*) "T ",n0*nc*dt," simulation time"
      

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
      alpha=w*mo*omega

c============================================================
c DIREKTER TERM 
c------------------------------------------------------------
c
c Summen auf Null
c
      sv=0.
      svv=0.
      q2av=0.
      xq2av=0.

c-----------------------------------------------------------
c Anfangspunkte
c
c 
c Oszillator
c
      do 25 la=1,a 
        do 23 j=1,3
             qi(j,la)=rgau()/(xo*sqrt(2.))
 23     continue
        do 24 j=1,3
             qf(j,la)=rgau()/(xo*sqrt(2.))
 24     continue
 25   continue 
c
c Projektil, Anfangswert: 
c Direkter Term: Pfad vom Zentrum1 zum Zentrum1
c
      do 20 j=1,3
             xi(j)=rgau()*b+qi(j,1)
20    continue
      do 21 j=1,3
             xf(j)=rgau()*b+qf(j,1)
21    continue
c--------------------------------------------------------------
c Anfangswert, straight line zwischen Anfangs und Endpunkt
c das kann man sicher besser machen !
c
c       print *,'Anfangswert fuer Langevin Gleichung machen :'
        do 100 i=0,(nb+1)
           do 110 j=1,3
              x(j,i)=(xf(j)-xi(j))*i/real(nb+1)+xi(j)
	      do 111 la=1,a
                 q(j,i,la)=(qf(j,la)-qi(j,la))*i/real(nb+1)+qi(j,la)
111           continue
110        continue
100     continue



c-------------------------------------------------------------
c Langevin Simulation
c------------------------------------------------------------- 
c Gleichgewichtslauf
c
c         print *,'Gleichgewichtslauf :',ia
          do 190 l2=1,ng 
c-----------------------
c ein schritt
c
             call makran(dw,dwq,nb,ds,a)
             do 191 i=0,nb+1
                call sd(i,x,q,z,zq,nb,r2,r2t)
                do 192 j=1,3
                   dx(j,i)=z(j)*dt+dw(j,i)
		   do 195 la=1,a
                      dq(j,i,la)=zq(j,la)*dt+dwq(j,i,la)
195                continue
192             continue
191          continue
             do 193 i=0,nb+1
                do 194 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
		  do 196 la=1,a
                      q(j,i,la)=q(j,i,la)+dq(j,i,la)
196               continue
194             continue
193          continue   
c
c ende des schrittes
c------------------------
190     continue
c----------------------------------------------------------
c eigentliche simulation (hier wird gemessen!!)
c
c       print *,'Simulationsbeginn ',ia
        do 200 l1=1,nc
c
c auslassung (n0 schritte)
c
          do 210 l2=1,n0 
c-------------------------------
c ein schritt
c 
             call makran(dw,dwq,nb,ds,a)
             do 220 i=0,nb+1
                call sd(i,x,q,z,zq,nb,r2,r2t)
                do 230 j=1,3
                  dx(j,i)=z(j)*dt+dw(j,i)
		  do 221 la=1,a
                     dq(j,i,la)=zq(j,la)*dt+dwq(j,i,la)
221               continue
230             continue
220          continue
             do 240 i=0,nb+1
                do 250 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
		  do 251 la=1,a
                     q(j,i,la)=q(j,i,la)+dq(j,i,la)
251               continue
250             continue
240          continue   
c
c ende des schrittes
c------------------------------
210     continue
c
c Potentialanteil der Wirkung berechnen
c
        st=0.
        do 400 i=1,nb
	    do 401 la=1,a
               st=st+v(r2(i,la))*db
401         continue
            st=st-vf(r2t(i))*db
400     continue
        do 402 la=1,a	
	   st=st+v(r2(nb+1,la))*db/2
	   st=st+v(r2(0,la))*db/2
402     continue
        st=st-vf(r2t(0))*db/2
        st=st-vf(r2t(nb+1))*db/2

        sv=sv+exp(st) 
	svv=svv+exp(2*st)
c
c mittlerer Abstand zwischen zwei Targets und zum Projektil am Halbpfad
c
	if (a.gt.1) then
	  do 410 i=1,3
		q2av=(q(i,nb/2,1)-q(i,nb/2,2))**2+q2av
		xq2av=(x(i,nb/2)-q(i,nb/2,1))**2+xq2av
410       continue
	endif

c
c evtl. protokollieren
c
	if ((mod(l1,10).eq.0).and.twrite) then
	   write(10,*) dt*l1, sv/l1, q2av/l1, xq2av/l1
        endif
c
c Konfiguration auf file schreiben
c
        if (cwrite) then
          do 300 i=0,nb+1
              write(10,8010) (x(j,i),j=1,3)
300       continue
        endif
c
c Pfadmittelpunkte schreiben
c
	if ((mod(l1,100).eq.0).and.hwrite) then
	   do 301 i=1,a
		write(15,8011) dt*l1,i,(q(j,nb/2,i),j=1,3)
301	   continue
	endif

	
200    continue
c
c hier endet die Simulationsschleife ueber nc Konfigurationen
c
c berechne mittleren potentialanteil der Wirkung

       sv=sv/nc
       svv=(svv/nc-sv**2)
       q2av=q2av/nc
       xq2av=xq2av/nc

       sad=sv
       savd=svv

c============================================================
c INDIREKTER TERM 
c------------------------------------------------------------
c Schleife ueber Anfangskonfigurationen
c
c Summen auf Null
c

      if (a.eq.1) goto 5000

      sv=0.
      svv=0.
      q2avi=0.
      xq2avi=0.

c-----------------------------------------------------------
c Anfangspunkte
c
c     print *,'Anfangswerte samplen :'
c 
c Oszillator
c
      do 1025 la=1,a 
        do 1023 j=1,3
             qi(j,la)=rgau()*xo
 1023   continue
        do 1024 j=1,3
             qf(j,la)=rgau()*xo
 1024     continue
 1025   continue 
c
c Projektil, Anfangswert: 
c Direkter Term: Pfad vom Zentrum1 zum Zentrum1
c
      do 1020 j=1,3
             xi(j)=rgau()*b+qi(j,1)
1020  continue
      do 1021 j=1,3
             xf(j)=rgau()*b+qf(j,1)
1021  continue
c--------------------------------------------------------------
c Anfangswert, straight line zwischen Anfangs und Endpunkt
c das kann man sicher besser machen !
c
c       print *,'Anfangswert fuer Langevin Gleichung machen :'
        do 1100 i=0,(nb+1)
           do 1110 j=1,3
              x(j,i)=(xf(j)-xi(j))*i/real(nb+1)+xi(j)
	      do 1111 la=1,a
                 q(j,i,la)=(qf(j,la)-qi(j,la))*i/real(nb+1)+qi(j,la)
1111          continue
1110       continue
1100    continue
c-------------------------------------------------------------
c Langevin Simulation
c------------------------------------------------------------- 
c Gleichgewichtslauf
c
c         print *,'Gleichgewichtslauf :',ia
          do 1190 l2=1,ng 
c-----------------------
c ein schritt
c
             call makran(dw,dwq,nb,ds,a)
             do 1191 i=0,nb+1
                call sid(i,x,q,z,zq,nb,r2,r2t)
                do 1192 j=1,3
                   dx(j,i)=z(j)*dt+dw(j,i)
		   do 1195 la=1,a
                      dq(j,i,la)=zq(j,la)*dt+dwq(j,i,la)
1195               continue
1192            continue
1191         continue
             do 1193 i=0,nb+1
                do 1194 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
		  do 1196 la=1,a
                      q(j,i,la)=q(j,i,la)+dq(j,i,la)
1196              continue
1194            continue
1193         continue   
c
c ende des schrittes
c------------------------
1190    continue
c----------------------------------------------------------
c eigentliche simulation (hier wird gemessen!!)
c
c       print *,'Simulationsbeginn ',ia
        do 1200 l1=1,nc
c
c auslassung (n0 schritte)
c
          do 1210 l2=1,n0 
c-------------------------------
c ein schritt
c 
             call makran(dw,dwq,nb,ds,a)
             do 1220 i=0,nb+1
                call sid(i,x,q,z,zq,nb,r2,r2t)
                do 1230 j=1,3
                  dx(j,i)=z(j)*dt+dw(j,i)
		  do 1221 la=1,a
                     dq(j,i,la)=zq(j,la)*dt+dwq(j,i,la)
1221              continue
1230            continue
1220         continue
             do 1240 i=0,nb+1
                do 1250 j=1,3
                  x(j,i)=x(j,i)+dx(j,i)
		  do 1251 la=1,a
                     q(j,i,la)=q(j,i,la)+dq(j,i,la)
1251              continue
1250            continue
1240         continue   
c
c ende des schrittes
c------------------------------
1210    continue
c
c Potentialanteil der Wirkung berechnen
c
        st=0.
        do 1400 i=1,nb
	    do 1401 la=1,a
               st=st+v(r2(i,la))*db
1401        continue
            st=st-vf(r2t(i))*db
1400    continue
        do 1402 la=1,a	
	   st=st+v(r2(nb+1,la))*db/2
	   st=st+v(r2(0,la))*db/2
1402    continue
        st=st-vf(r2t(0))*db/2
        st=st-vf(r2t(nb+1))*db/2

        sv=sv+exp(st) 
	svv=svv+exp(2*st)

c
c mittlerer Abstand zwischen zwei Targets und zum Projektil am Halbpfad
c
	if (a.gt.1) then
	  do 1410 i=1,3
		q2avi=(q(i,nb/2,1)-q(i,nb/2,2))**2+q2avi
		xq2avi=(x(i,nb/2)-q(i,nb/2,1))**2+xq2avi
1410       continue
	endif

c
c evtl. protokollieren
c
	if ((mod(l1,10).eq.0).and.twrite) then
	   write(11,*) dt*l1,sv/l1, q2avi/l1, xq2avi/l1
        endif
c
c Konfiguration auf file schreiben
c
        if (cwrite) then
          do 1300 i=0,nb+1
              write(10,8010) (x(j,i),j=1,3)
1300      continue
        endif
1200    continue
c
c hier endet die Simulationsschleife ueber nc Konfigurationen
c
c berechne mittleren potentialanteil der Wirkung
c
       sv=sv/nc
       svv=(svv/nc-sv**2)
       xq2avi=xq2avi/nc
       q2avi=q2avi/nc
c 
c Der indirekte Term
c

       said=sv
       savid=svv

5000   continue
c---------------------------------------------------------
c ENDE INDIREKTER TERM
c---------------------------------------------------------
c Berechnung der Querschnittsverhaeltnisse
c----------------------------------------------------------
       
       if (a.gt.1) then
       	   sr=1./a*1./sad+(a-1.)/a*1./said
       else
       	   sr=1./a*1./sad
       end if 

       write (9, *) "sigmaratio ", sr, " ratio of cross sections"
       write (9, *) "direct ", sad, " inverse ratio from direct"
       write (9, *) "errdirect ", savd, " variance of direct"
       write (9, *) "xq2av ", xq2av,  " projectile target mean square"
       write (9, *) "q2av ", q2av,  " target target mean square"
       if (a.gt.1) then
          write (9, *) "indirect ", said, " inverse ratio from indirect"
          write (9, *) "errindirect ", savid, " variance of indirect"
          write (9, *) "xq2avi ", xq2avi,  " projectile target indirect"
          write (9, *) "q2avi ", q2avi,  " target target indirect"
       end if
       write (9, *) "</RESULT>"

       write(14,*) (nb+1)*db,1/((nb+1)*db),sad,said,sr
       write(16,*) (nb+1)*db,1/((nb+1)*db),q2av,xq2av
       write(12,*) (nb+1)*db,1/((nb+1)*db),savd,sqrt(savd),
     A             savid,sqrt(savid)


       write(9,*) "</SIMULATION>"
       goto 1


8010   format(3e12.4)
8011   format(e12.4,i8,3e12.4)
9999   continue
       end
      
      subroutine sd(i,x,q,z,zq,nb,r2,r2t)
c
c berechnet die ableitung der Wirkung am iten Schritt
c
      implicit none
      integer nbmax,amax
      parameter(nbmax=8000,amax=20)
      real b,v0,db,m,pi,w,mb,mob,f1,f2
      integer i,j,nb,a,la
      real r2(0:nbmax+1,amax),r2t(0:nbmax+1)
      real x(3,0:nbmax+1),q(3,0:nbmax+1,amax)
      real z(3),zq(3,amax)
      real vp(amax),r(amax),rt
      real omega,mo,xo
      real alpha
      common/const/pi
      common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2,a,alpha
 
      do 1 la=1,a
         r(la)=0. 
1     continue
      rt=0.
      do 10 j=1,3
	 do 11 la=1,a
            r(la)=r(la)+(x(j,i)-q(j,i,la))**2
11       continue
	 rt=rt+x(j,i)**2
10    continue
      r2t(i)=rt
      do 12 la=1,a
         r2(i,la)=r(la)
         vp(la) =-2./w*v0*exp(-r(la)/w)*db
12    continue
      if (i.eq.0) then
        do 100 j=1,3
c Projectile:
          z(j)=mb*(x(j,i+1)-x(j,i))
     A          -(x(j,i)-q(j,i,1))/b**2
	  do 101 la=1,a
              z(j)=z(j)-vp(la)*(x(j,i)-q(j,i,la))/2.
101       continue
c Target:
	  do 102 la=1,a
             zq(j,la)=f1*(q(j,i+1,la)-f2*q(j,i,la))
     A         -vp(la)*(q(j,i,la)-x(j,i))/2.
     B         -xo**2*q(j,i,la)
102       continue
	  zq(j,1)=zq(j,1)-(q(j,i,1)-x(j,i))/b**2
100     continue
      elseif (i.eq.(nb+1)) then
        do 110 j=1,3
c Projectile: 
          z(j)=mb*(-x(j,i)+x(j,i-1))
     A        -(x(j,i)-q(j,i,1))/b**2
	  do 112 la=1,a
	      z(j)=z(j)-vp(la)*(x(j,i)-q(j,i,la))/2.
112       continue
c Target: 
	  do 111 la=1,a
             zq(j,la)=f1*(-f2*q(j,i,la)+q(j,i-1,la))
     A            -vp(la)*(q(j,i,la)-x(j,i))/2.
     B            -xo**2*q(j,i,la)
111       continue
	  zq(j,1)=zq(j,1)-(q(j,i,1)-x(j,i))/b**2
110     continue
      else
        do 120 j=1,3
c Projectile:     
          z(j)=mb*(x(j,i+1)-2.*x(j,i)+x(j,i-1))
	  do 122 la=1,a
	      z(j)=z(j)-vp(la)*(x(j,i)-q(j,i,la))
122       continue
c Target:
	  do 121 la=1,a
             zq(j,la)=f1*(q(j,i+1,la)-2.*f2*q(j,i,la)+q(j,i-1,la))
     B            -vp(la)*(q(j,i,la)-x(j,i))
121       continue
120     continue
      endif
      end

      subroutine sid(i,x,q,z,zq,nb,r2,r2t)
c
c berechnet die ableitung der Wirkung am iten Schritt
c
      implicit none
      integer nbmax,amax
      parameter(nbmax=8000,amax=20)
      real b,v0,db,m,pi,w,mb,mob,f1,f2
      integer i,j,nb,a,la
      real r2(0:nbmax+1,amax),r2t(0:nbmax+1)
      real x(3,0:nbmax+1),q(3,0:nbmax+1,amax)
      real z(3),zq(3,amax)
      real vp(amax),r(amax),rt
      real omega,mo,xo
      real alpha
      common/const/pi
      common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2,a,alpha
 
      do 1 la=1,a
         r(la)=0. 
1     continue
      rt=0.
      do 10 j=1,3
	 do 11 la=1,a
            r(la)=r(la)+(x(j,i)-q(j,i,la))**2
11       continue
	 rt=rt+x(j,i)**2
10    continue
      r2t(i)=rt
      do 12 la=1,a
         r2(i,la)=r(la)
         vp(la) =-2./w*v0*exp(-r(la)/w)*db
12    continue
      if (i.eq.0) then
        do 100 j=1,3
c Projectile:
          z(j)=mb*(x(j,i+1)-x(j,i))
     A          -(x(j,i)-q(j,i,1))/b**2
	  do 101 la=1,a
              z(j)=z(j)-vp(la)*(x(j,i)-q(j,i,la))/2.
101       continue
c Target:
	  do 102 la=1,a
             zq(j,la)=f1*(q(j,i+1,la)-f2*q(j,i,la))
     A         -vp(la)*(q(j,i,la)-x(j,i))/2.
     B         -xo**2*q(j,i,la)
102       continue
	  zq(j,1)=zq(j,1)-(q(j,i,1)-x(j,i))/b**2
100     continue
      elseif (i.eq.(nb+1)) then
        do 110 j=1,3
c Projectile: 
          z(j)=mb*(-x(j,i)+x(j,i-1))
     A        -(x(j,i)-q(j,i,2))/b**2
	  do 112 la=1,a
	      z(j)=z(j)-vp(la)*(x(j,i)-q(j,i,la))/2.
112       continue
c Target: 
	  do 111 la=1,a
             zq(j,la)=f1*(-f2*q(j,i,la)+q(j,i-1,la))
     A            -vp(la)*(q(j,i,la)-x(j,i))/2.
     B            -xo**2*q(j,i,la)
111       continue
	  zq(j,2)=zq(j,2)-(q(j,i,2)-x(j,i))/b**2
110     continue
      else
        do 120 j=1,3
c Projectile:     
          z(j)=mb*(x(j,i+1)-2.*x(j,i)+x(j,i-1))
	  do 122 la=1,a
	      z(j)=z(j)-vp(la)*(x(j,i)-q(j,i,la))
122       continue
c Target:
	  do 121 la=1,a
             zq(j,la)=f1*(q(j,i+1,la)-2.*f2*q(j,i,la)+q(j,i-1,la))
     B            -vp(la)*(q(j,i,la)-x(j,i))
121       continue
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
      real alpha
      integer a
      common/const/pi
      common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2,a,alpha

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
      integer a
      common/const/pi
      common/par/m,b,v0,db,w,omega,mo,xo,mb,mob,f1,f2,a,alpha

      vf=v0*a*sqrt(alpha/(1+alpha))**3*exp(-r/w*alpha/(1+alpha))

      end


      
      REAL FUNCTION RAN()

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE CALCULATES A RANDOM NUMBER.                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      DOUBLE PRECISION X,XK,PM,A,C,IM
      LOGICAL BUILDIN
      PARAMETER(A=14662125.D00,C=13136923.D00,PM=2.D00**48)
      PARAMETER(IM=2147483647)

      COMMON/RG/X,XK,BUILDIN

      if (buildin) then
	ran=real(irand(0))/IM
      else
        XK=XK+.5D00
        X=A*X+DINT(C*XK)
        X=X-DINT(X/PM)*PM
        RAN=REAL(X/PM)
      endif
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


      subroutine makran(w,wq,n,ds,a)
      implicit none
c
c maximale zahl der beta schritte
c
      integer nbmax,amax
      parameter(nbmax=8000,amax=20)
c
c felder fuer zufallszahlen
c
      real w(3,0:nbmax+1)
      real wq(3,0:nbmax+1,amax)
      real x1,x2
      real cs,sn,lg
      real ds
      real ran
      external ran
      integer n,a,la
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
	   do 110 la=1,a
111          x1=ran()
	     if (x1.lt.1e-16) goto 111
             x2=ran()
             cs=cos(2.*pi*x2)
             sn=sin(2.*pi*x2)
             lg=sqrt(-2.*alog(x1))*ds
             wq(j,i,la)  =lg*cs
             wq(j,i+1,la)=lg*sn
110   continue
      end
