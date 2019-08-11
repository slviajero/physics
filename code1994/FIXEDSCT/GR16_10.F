      PROGRAM green16 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  ORGANISATION OF THE PROGRAM:                                        C
C    SUBROUTINE EDCS  CALCULATES GREENS FUNCTION FOR A GIVE CONFIGURATION 
C                      FUNCTIONS                                       C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   PARAMETER USED IN THE WHOLE PROGRAM:                               C                          C
C      M      DIMENSIONALITY OF THE SPACE                              C
C      N      NUMBER OF PARTICLES                                      C
C      B      WIDTH OF THE HARMONIC OSCILLATOR WAVE FUNCTION           C
C   VARIABLES USED IN THE STOCHASTIC PART OF THE PROGRAM               C
C      NSTART  NUMBER OF INITIAL CONFIGURATIONS                        C
C      NKEEP   NUMBER OF CONFIGURATION KEPT AFTER EACH INITIAL         C
C              CONFIGURATION                                           C
C      X       INITIAL VALUE FOR RANDOM NUMBERS                        C
C      XD      COUNTER FOR RANDOM NUMBERS                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c
c Version string
c
      character*(80) version
      
      
      COMMON/PARA/PI,NSTART,NKEEP

      PI=4.*ATAN(1.)
      version='Version 1.01 7.5.1992       '
      NAUS=1

      OPEN (1,FILE='data',STATUS='UNKNOWN')
      OPEN (2,FILE='curve',STATUS='UNKNOWN')
      OPEN (3,FILE='CONFIG',STATUS='UNKNOWN')
      open (4,FILE='report')
      open (16,file='f1data')

      REWIND(1)
      REWIND(2)
      REWIND(3)
      rewind(4)

      write(4,*) 'Programm Green:'
      write(4,*) version
      write(4,*) 'Plot1 : <Re G(t)> '
      write(4,*) 'Plot2 : <Im G(t)> '
      write(4,*) 'Plot3 : (Abs <G(t)>)**2 '
      write(4,*) 'Plot4 : < (Abs G(t))**2 > '
      write(4,*) 'Plot5 : (Abs G0(t))**2    '
      write(4,*) 'Plot6 : Im < G(E) >       '
      write(4,*) 'Plot7 : Im < G0(E) >      '

      READ(3,1010,ERR=100,END=100)NSTART,NKEEP
        CALL EDCS_s
300   CLOSE (1,STATUS='KEEP')
      CLOSE (2,STATUS='KEEP')
      CLOSE (3,STATUS='KEEP')

100   continue
1010  format(2i8)      
      END


      SUBROUTINE EDCS_s
c
c
c   calculates the greensfunction g(t,x,y) at x=y=0
C   first, the energy dependent greensfunction g(e,x=y=0)-g0(e,x=y=0)
c   is calculated for a number of energies (the free greens function
c   is singular an dtherefore subtracted)
c   then the result is fourier transformed and writem to a file
c   for the calculation of the energy dependent greensfunction
c   see D. Stoll
c----------------------------------------------------------------------
c
c Program history:
c
c 7.5.1992: old program without time averaging is
c           tested after system crash, averaging 
c           of energy dependent greensfunction is   
c           included
c 11.5.1992: input of scattering amplitudes modifies
c            phaseshift can be given
c 12.5.1992: abs of averaged greensfunction is
c            given
c            logarithimic timescale included
c
C   PARAMETER:                                                         C
C      KAMAX   NUMBER OF VALUES OF THE SCATTERING ANGLE WHERE THE      C
C              DIFFERENTIAL CROSS SECTION IS CALCULATED                C
C      KMAX    UPPER LIMIT FOR THE NUMBER OF VALUES OF THE INCIDENT    C
C              MOMENTUM                                                C
C      Q0      UPPER LIMIT FOR CALCULATION OF SCATTERING AMPLITUDE IN  C
C              FORWARD DIRECTION ONLY                                  C
C   IMPORTANT VARIABLES:                                               C
C      AM      CONTAINS THE MATRIX ELEMENTS OF THE MULTIPLE SCATTERING C
C              EQUATIONS                                               C
C      NTHET   ANGULAR STEP                                            C
C      QFACT   INCIDENT MOMENTUM STEP                                  C
C      FKKP    CONTAINS THE SCATTERING AMPLITUDES                      C
C      MULTI   MULTIPLICITY OF A CONFIGURATION                         C
C      PHI     INHOMOGENEOUS VECTOR OF THE MULTIPLE SCATTERING         C
C              EQUATIONS IN THE BEGINNING SOLUTION OF THESE EQUATIONS  C
C              IN THE END                                              C
C      PHI1    CONTAINS THE OUTGOING WAVE                              C
C      Q       MOMENTUM IN Z-DIRECTION                                 C
C      QANGLE  DIRECTIONS OF THE OUTGOING WAVE                         C
C      RMS     RMS-RADIUS OF THE DENSITY                               C
C      SLENG   PROJECTILE NUCLEON SCATTERING LENGTH                    C
C      SUM     TOTAL NUMBER OF CONFIGURATIONS                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER(M=3,N=16,KAMAX=61,KMAX=10,Q0=.2,kfm=10,ne=200,ntst=50)
      parameter(nint=400)
c
c   variables which are also used in the scattering program
c

      REAL FIELD(M,N),RIJ(N,N)
      COMPLEX IC,AM(N,N),PHI(N),FAKT,
     1        PHI1(N),sleng0,pleng(kmax),sleng(kmax)

c
c    new ones
c
      real ri(n), energ(0:ne)
c
c  non averaged greensfunction variables 
c
      complex*16 greene(kmax,0:ne),greent(kmax,ntst)
c
c  averaged greensfuntion variables
c
      complex*16 grav(kmax,ntst)
      real grabs(kmax,ntst),grave(kmax,0:ne)      
      real tv(ntst),tres
c
c this one is for the fourier integration of the greensfunction
c it contains the linearly interpolated energy greensfunction
c
      complex*16 interp(nint),interu(nint)
      complex*16 dg,gn,dgdk,gsum
      real dk
c
c this is the phase shift flag, if it is true, then the
c second two lines of f1data contain the phaseshift
c of the s-wave scattering amlitude. The scattering amplitude
c is then generated as
c f = ( exp 2i delta - 1 ) / 2ik
c
c
      logical phase

      COMMON/PARA/PI,NSTART,NKEEP

c
c important parameters
c
       tres  = 0.015
       btres = 1.2
c
C  INPUT OF THE ANGULAR STEPS, MOMENTUM AND P-WAVE SCATTERING LENGTH
c
c  this program uses the same input an output files as
c  the scattering programs 
c  the momentum q is now the maximum energy for
C the fourier transform

      READ(16,1045)KNUMB
      read(16,*)sleng
      READ(16,*)PLENG
      NTHET=180
      READ(16,*)Q
      READ(16,*)NTHET
      NTHET=MAX0(4,NTHET)
      mthet=nthet
      kangle=181
      read(16,*)kf
c
c  kf is used to determine phase
c
      phase=(kf.eq.1) 


c
c     here the energy values for the fourier transform 
c     are determined.
c
c

      qm=q
      dk=qm/ne
      do 7 i=0,ne
	energ(i)=qm/ne*i
7     continue
c
c     here the time values for the fourier transform 
c     are determined
c     this is not done very carefully, of course the
c     resolution has to be compared with the energy 
c     resolution
c     here, a logarithmic timescale is used
c     to extend the range of time values
c
      do 8 i=1,ntst
	 tv(i)=tres*(btres)**i
8     continue

      IC=CMPLX(0.,1.)
      SUM=0.
      RMS=0.
c
c load configurations (nstart*nkeep)
c store one configuration in FIELD
c and solve linear equations
c
      nc=nstart*nkeep
      DO 200 K1=1,NSTART*NKEEP
        READ(3,*,END=800)FIELD
        READ(3,1010,END=800)RSQU,MULTI
c
C  BEGINNING OF THE CALCULATION OF THE SCATTERING AMPLITUDE FOR ONE
C  CONFIGURATION
        SUM=SUM+MULTI
        RMS=RMS+RSQU*MULTI
c
C  DETERMINATION OF THE PARTICLE DISTANCES
c  absolute of r_i - r_j is calculated
c  and stored in rij
c
	DO 20 K2=1,N-1
          DO 20 K3=K2+1,N
            R=0.
            DO 10 K4=1,M
              R=R+(FIELD(K4,K2)-FIELD(K4,K3))**2
10          CONTINUE
            RIJ(K2,K3)=R**.5
20      CONTINUE
c
c  absolute of r_i is calculated
c  and stored in ri
c
	DO 21 K2=1,N
            R=0.
            DO 11 K4=1,M
              R=R+(FIELD(K4,K2))**2
11          CONTINUE
            RI(K2)=R**.5
21      CONTINUE
c
c start of calculation of the greens function
c first, the inhomogenous vector is calculated,
c then the system of linear equations is
c solved.

c
c this loop controls the different scattering amplitude
c
       DO 200 K2=1,KNUMB
c
c this loop controls the different energies
c
          do 120 nei=0,ne
               q=energ(nei)
c 
c according to the value of phase, the
c scattering amplitude is generated
c as a unitary scattering amplitude
c or simply taken from sleng
c
	    if (phase) then
	      if (q.ne.0) then
	       sleng0=(exp(2*ic*pleng(k2))-1)/(2*ic*q)
              else 
	       sleng0=0
              endif
            else
              sleng0=sleng(k2)
            endif
c
C   CALCULATION OF THE INHOMOGENEOUS VECTORS
c   for the greens function calculation this is just
c   h_0(k*ri) the scattering amplitude is multiplied to 
c   the outgoing vector
c   keep in mind, that this algorithm is very stupid,
c   because the same calculation has to be done 
c   for the different scattering amplitudes, but
c   this is fortunatly not very time consuming 
c
          DO 30 K3=1,N
            PHI (K3)  =CEXP(IC*Q*ri(K3))/ri(k3)
            PHI1(K3)  =PHI(K3)*sleng0
30      CONTINUE

c
c
c here, the greensfunction matrix is build up
c
c
C   CALCULATION OF THE MATRIX ELEMENTS
          DO 40 K3=1,N-1
            AM(K3,K3)=1.
            DO 40 K4=K3+1,N
              AM(K3,K4)=-SLENG0*CEXP(IC*Q*RIJ(K3,K4))/RIJ(K3,K4)
              AM(K4,K3)=AM(K3,K4)
40        CONTINUE
          AM(N,N)=CMPLX(1.)
c
c
c  some test aids , not used any longer
c
c


c         do 7777 iii=1,n
c           do 7777 jjj=1,n
c777      write(1,*) am(iii,jjj)


c         write(1,*)'Confi'
c         write(1,*)field
c         write(1,*)'matrix'
c         write(1,*)am
c         write(1,*)'inhomo'
c         write(1,*)phi


c
c linear equations are solved here with an algorithm, which uses
c the symmetry of the matrix
c

C   SOLUTION OF THE SYSTEM OF MULTIPLE SCATTERING EQUATIONS
          DO 70 K3=1,N-1
            DO 70 K4=K3+1,N
              FAKT=AM(K3,K4)/AM(K3,K3)
              AM(K3,K4)=FAKT
              DO 70 K5=K3+1,N
                AM(K5,K4)=AM(K5,K4)-AM(K5,K3)*FAKT
70        CONTINUE
          DO 90 K4=2,N
            DO 90 K5=1,K4-1
              PHI(K4)=PHI(K4)-AM(K5,K4)*PHI(K5)
90        CONTINUE
          PHI(N)=PHI(N)/AM(N,N)
          DO 110 K4=1,N-1
            K3=N-K4
            DO 100 K5=K3+1,N
              PHI(K3)=PHI(K3)-AM(K5,K3)*PHI(K5)
100         CONTINUE
            PHI(K3)=PHI(K3)/AM(K3,K3)
110       CONTINUE

c         write(1,*)'loesung'
c         write(1,*)phi

c
c phi contains now the solution for inhomogenous vector
c corresponding to scattering length k2 and energy nei
c now this solution is projected on phi1 to give the energy dependent
c greens function minus the free one
c

          greene(k2,nei)=(0.,0.)
          do 120 k6=1,n
	     greene(k2,nei)=greene(k2,nei)+phi1(k6)*phi(k6)
120       continue
	  
c
c also the loop which counts the energies ends here.    
c here, the variable greene contains the energy dependent
c greens function for one configuration and one scattering length
c

c now the time dependent greensfunction is calculated
c all this assumes that there are no E<0 singularities
c such that only positive energies contribute to the fourier
c integral
c
c note that energ is the absolute of the wavevector, which
c is the square of the energy. thus in the fourier integral, the
c density of states has to be taken into  account explicitly
c
c this is done via the identity
c int dE G(E) exp (-iEt) = int dk 2k G(k) exp( -i k*kt)
c
c 
c for all time values
c         do 130 i=1,ntst
c            greent(k2,i)=(0.,0.)
c            do 130 j=0,ne
c             q=energ(j)
c             greent(k2,i)=greent(k2,i)+(qm/ne)*2*q*cexp(-ic*q*q*tv(i))*
c    A                                    greene(k2,j)
c130       continue
c
c
c the above method of doing the fourier transform
c is not to clever, because of the strong oszillations
c of exp -i k*k t compared to greene
c
c the energy greensfunction is linearly interploated
c first up to nint points
c this save computation time 
c note, that the last point has to be done by hand
c
          do 131 i=0,ne-1
	    dg=greene(k2,i+1)-greene(k2,i)
	    gn=greene(k2,i)
	    dgdk=dg/dk
	    do 131 j=0,(nint/ne)-1
              ie=(nint/ne)*i+j
	      interu(ie)=dgdk*(qm/nint*ie-energ(i))+gn
131       continue
	  interu(nint)=greene(k2,ne)
c
c after the linear interploation all values are multiplied
c with the exponential and the integration measure
c and then a simpsom integration is performed
c 
c
c for all times
c
	 do 134 jjj=1,ntst,1
c
c calculate integrand
c
	  do 132 i=0,nint
             q=qm/nint*i
	     interp(i)=(qm/nint)*2*q*cexp(-ic*q*q*tv(jjj))*interu(i)
132       continue
c
c integrate
c
          gsum=(0.,0.)
          do 133 i=1,nint-1,2
	    gsum=gsum+4*interp(i)+2*interp(i+1)
133       continue
          gsum=(gsum+interp(0)-interp(nint))/3
c
c store result
c
         greent(k2,jjj)=gsum
134      continue
c
c
c now greent contains the time dependent greens function
c for one configuration and one value of sleng
c
c at this time, averaging over configurations has to de done
c we average the greensfunction and its absolute squared
c but first, the free greensfunction has to be added
c for factors, see Economu
c
c
c
c first the time dependent greensfunction is
c averaged
c
       do 140 i=1,ntst
         grav(k2,i)=grav(k2,i)+greent(k2,i)
     A               +1/(sqrt(ic*4*pi*tv(i))*tv(i))
          grabs(k2,i)=grabs(k2,i)
     A      +abs(      
     B        greent(k2,i)
     C       +1/(sqrt(ic*4*pi*tv(i))*tv(i))
     D       )**2
140    continue
c
c here the energy dependent greensfunction is averaged
c as we calculate g(E,x,x) only the imaginary part
c of the free greensfunction is well defined, so
c only the imaginary part is stored in grave 
c
       do i=0,ne-1
         grave(k2,i)=grave(k2,i)+imag(greene(k2,i))
     A               +sqrt(energ(i))
       end do


c
c now the calculation is done for the next scattering amplitude
c and then for the next configuration
c
200   CONTINUE
c
c
c now the information is written to a file
c yhis is not done very carfully at this stage     
c
      do 300 i=1,knumb
	if (phase) then
	   write(4,*) 'Phaseshift : ',pleng(i)
        else
           write(4,*) 'Scattering amplitude :',sleng(i)
        endif
        write(4,*) 'Plotnummer : ',i
	do j=1,ntst
	 write(2,1100) 1+i*10,tv(j),real(grav(i,j))/nc
	 write(2,1100) 2+i*10,tv(j),imag(grav(i,j))/nc
	 write(2,1100) 3+i*10,tv(j),(abs(grav(i,j))/nc)**2
	 write(2,1100) 4+i*10,tv(j),grabs(i,j)/nc
	 write(2,1100) 5+i*10,tv(j),abs(1/(sqrt(ic*4*pi*tv(j))*tv(j)))**2
        end do
	do j=0,ne-1
	 write(2,1100) 6+i*10,energ(j),grave(i,j)/nc
	 write(2,1100) 7+i*10,energ(j),sqrt(energ(j))
        end do
300   continue

800   RETURN


1005  FORMAT(E12.4)
1010  FORMAT(E12.4,I8)
1045  FORMAT(I3)
1100  FORMAT(I3,2E12.4)
1110  FORMAT(I3,4E12.4)

      END
