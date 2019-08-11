	Program Osc10
C	Zufallszahlengenerator korrigiert
C	Osz. alleine, exakter Kurzzeitpropagator
C	Ergebnis: 3/MT(Fr+OM), Norm=1!!
C	Mit Verdopplung der Konfiguration bei starkem Schwund
C	============
C	============
        implicit none

C	Felder
C	======

        Integer Parco1, Parco2
	Parameter(Parco2=70000)
	Parameter(Parco1=100000)
	Real XT(Parco1,3), YT(Parco1,3)
	Real YYT(Parco1,3)
	Real XT2, YT2, DX2, DY2, DY3, DY4

C	Zaehler
C	=======

	Integer  Step, Nstep
	Integer  Conf, Nconf, D
	Integer Count,Count1,Count2,Anf1,Anf2,Cou

C	Abfrage/Sonstige
C	================
	
	Real W, DW, WN, WZ,H
	Integer Wint,Wint2

C	Zufallszahlen
C	=============

	Integer K1,ZufC
	Real PI, TT,ZUFZ
	Real Ran, Gau
	External Ran, Gau
	Double Precision XR, XK
	Common/RG/XR,XK
	Common/Const/Pi

C	Eingabeparameter
C	================

	Real MT, V ,Fr ,Om,Om2, Ommt, Ommt2, eps,rr,r

C	Normierungen
C	============

	Real NN,Nosc,NT1,NT2,N0, E0,DifT

C	Initialisierung
C	===============
C	Pi
	Pi=4.*atan(1.)
	open(10,file='dat4')


C	Zufallszahlengenerator
C	======================

	XR=PI*1.D11
	XK=0.D00
	DO 10 K1=1,100000
	TT=TT+Ran()
10	Continue
	print *,TT/100000

C	Eingabe
C	=======

	print *,'Zahl der Schritte'
	read *,Nstep

	print *,'Schrittweite'
	read *,Eps


	print *,'Masse Target'
	read *,MT

	print *,'Frequenz'
	read *,Om,Fr
	

C	Normierungen
C	============
	
	NT1=sqrt(1.0/(Om*MT))
C	Nosc=2.*Pi*sqrt(MB)*V0*V0*2.**(3./2.)*(b**3)

	DifT=sqrt(sinh(eps*OM)/MT/OM)

        E0=(3./2.)*Om


        Ommt=Om*MT	
	Ommt2=Ommt**2
	Om2=OM*OM

C	Anfangskonfigurationen
C	======================

	Nconf=Parco2

	do 50 Conf=1,Nconf
	YT2=0
	XT2=0
		do 40 D=1,3
		 XT(Conf,D)=1./sqrt(Fr*MT)*Gau()
40		continue

50	continue


C	Propagation/Auswahl
C	===================


	do 1000 Step=1,Nstep
           
	      Anf1=0
              Anf2=0

	      do 610 Conf=1,Nconf

                   XT2=0

	           do 600 D=1,3
		      XT2=XT2+XT(Conf,D)**2
600	           continue

                      V=MT*Om*(XT2)*(cosh(eps*OM)-1.)/sinh(eps*OM)
                      
                      V=(V-eps*E0)/2.

		      if (V.ge.50) then
		       W=0.
                      else
		       W=exp(-V)
                      endif 
		    
		      Wint=int(W)
		      DW=(W-Wint)*.9994

		      if (Ran().lt.DW) then
			 Wint=Wint+1
                      endif
        

                      if (Wint.gt.0) then

	                do 608 Count1=1,Wint

                           Anf1=Anf1+1
		           YT2=0
                          
                           do 607 D=1,3
			      YT(Anf1,D)=XT(Conf,D)+DifT*Gau()
			      YT2=YT2+YT(Anf1,D)**2
607			  continue	
       
                           
	                   V=MT*OM*YT2*((cosh(eps*OM)-1.)/sinh(eps*OM))
			   V=(V-eps*E0)/2.

			   if (V.ge.50) then
			    W=0.
			   else	
			    W=exp(-V)
                           endif 
           
			   Wint2=int(W)
			   DW=(W-Wint2)*.9994
			    

			   IF (Ran().lt.DW) then
			       Wint2=Wint2+1
                           Endif

			   If (Wint2.gt.0) then

			     Do 655 Count2=1,Wint2

			     Anf2=Anf2+1

                               Do 653 D=1,3
				YYT(Anf2,D)=YT(Anf1,D)
653			       continue
655			     continue

			   Endif
                     


608                     continue

	              endif
C		Abfrage Wint nicht 0 endet hier!!!


610	       continue 


	       Nconf=Anf2
	       Anf2=0
	       Anf1=0

               do 730 Conf=1,Nconf
		  do 720 D=1,3
		     XT(Conf,D)=YYT(Conf,D)
720		  continue
730	       continue

C	Ergebnis
C	========

	       W=0
	       WN=0
	       WZ=0

               print*,rr,Nconf 
	       if (mod(Step,10).eq.0) Then
		  do 760 Conf=1,Nconf

			XT2=0

			do 755 D=1,3
                          XT2=XT2+XT(Conf,D)**2
755	                continue

                        V=Fr*MT/2.*XT2

			if (V.ge.50) then
			  W=0.
                        else
			  W=exp(-V)
                        endif

			WZ=WZ+W*XT2*XT2
			WN=WN+W
760	       continue		

C	Ausgabe
C	=======
C
C	4. Moment
C	=========
	       	  W=WZ/WN
	

C		  write(10,*)MT,Fr,OM,eps,Parco2,Nconf
C                 write(10,*)step*eps,abs(w)
		  print *,abs(w),Nconf,Step*eps

C	 Gerechneter Wert
C	 ================

                  r=2.*eps*step*OM

		  if (r.ge.50) then
		   V=0.
                  else
		   V=exp(-r)
                  endif

	          W=(2.*OM*Fr)/(1.-V)
	          W=W/(OM2+FR*Fr+2.*Fr*Om*(1./tanh(step*eps*Om)))
	          W=W**(3./2.)

	          print*,W,'Calc'

C	Simulierter Wert
C	================

                  print*,WN/parco2,'Sim'
C                 WN=WN/parco2/cou*((FR+OM)**2./(2.*Fr*Om))**(3./2.)
                  write(10,*)W,eps*step	
                  write(10,*)Wn/parco2,eps*step

               endif

1000    continue   		      
	end

C
C       Unterprogramme
C	==============	

	real function Gau()

	real z1, z2,z3
	real Pi
	real ran
	external ran
	common/Const/pi
	
222	z1=ran()
	z2=ran()
	z3=10.**(-16)

	if ((z1.le.z3).or.(z2.le.z3)) goto 222
	
	 
	Gau=sqrt(-2.*alog(z1))*cos(2.*pi*z2)
	end
        
	real function Ran()
	real Ran1
	double precision X, XK, pm, a, c
	parameter(a=14662125.D00,c=13136923.D00,pm=2.D00**48)
        Common/RG/X,XK

	XK=XK+.5D00
	X=a*x+dint(c*XK)
	X=X-dint(X/pm)*pm
	Ran=real(X/pm)
	
	return
	end
