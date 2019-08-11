      program streu
C
C
C
C    Dieses Programm berechnet die Streuamplitude eines
C    Vielfachstreuproblems, wie es in
C
C            Boundary Condition Models for Multiple
C            Scattering off Composite Systems
C
C    beschrieben ist. Es kann zus tzlich ueber
C    Monte Carlo Integration eine Mittelung ueber
C    eine gegebene Verteilung von Streuzentren
C    berechnet werden.
C
C----------------------------------------------------------
C   Auxilliary variables of the main program
C
       integer n
       real dummy
       character*(80) la
       integer bstat
C
C----------------------------------------------------------
C   Common block for the random number generator
C
       double precision x,xk
       COMMON/RG/X,XK
C
C----------------------------------------------------------
C  Here the global variables of the program are defined.
C  An initial value is assigned.
C  If there is a file scatpar, the initial values are
C  taken  from it.
C
C  This is the maximum number of particles
C
       parameter(maxpnm=20)
C
C  All variables are defined ...
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi
C
C  And assign to common - blocks
C
       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C  Load file scatpar, if existent. If file empty,
C  use standard initial values (at 200)
C
       open(98,FILE='scatpar',STATUS='unknown')
       rewind(98)
C  File empty : goto 200
       read(98,*,END=200,ERR=200) pnum
          read(98,202) config
          read(98,202) prtfil
          read(98,*) makprt
          read(98,*) k
          read(98,*) f0
          read(98,*) f1
          read(98,*) width
          read(98,*) hardc
          read(98,*) nconf
          close(98)
 202      format(a20)
       goto 201
C
C  No file scatpar, standard initial values are used
C
 200   continue
C
C Number of particles set
       pnum = 2
C Configuration file name
       config  = 'CONFIG'
C Make protokol file 'prtfil' if makprt is on
       makprt  = .FALSE.
       prtfil  = 'PROTO'
C Energy of the projectile
       k  = 1.0
C Width of the distribution
       width = 1.0
C Hardcoreradius
       hardc = 0.0
C Number of configurations
       nconf = 100
C s+p-wave scattering amplitude of scatterer
       f0  =  (1.0 , 1.0 )
       f1  =  (0.0 , 0.0 )
C
C      Non loadable parameters
C
 201   continue
C File, which is used to store the results in
       result='RESULT'
C In this file(s), the results are written in
C a format, which is required by an ERLGRAPH
C plot program. The plots are numbered
C by a integr cplo
       plotf='PLO'
       cplo = 21
C
C  Graf file opened
C
       open(96,FILE='plo1')
       rewind(96)

C Number of angular steps ( maximum 180 = 1 degree resolution )
       thetas = 20
C Pi
       pi=3.14125
C Bath mode flag  ( if true, standard output is surpressed )
       batch=.FALSE.

C
C----------------------------------------------------------
C   INITIALIZATION OF THE RANDOM NUMBER GENERATOR
C
      X=3.1415D11
      XK=0.D00
      DO 99 N=1,20
        dummy=RAN()
99    CONTINUE
C
C----------------------------------------------------------
C  Open result file
C
       open(99,FILE=result)
C
C  In this mode, all new results are appended to the
C  old file. If C before REWIND is removed, a new
C  file at every run will be generated
C
C      rewind(99)
C
C----------------------------------------------------------
C  The program supports interactive or batch operation.
C  In this section the file scatbat is opened. If the file
C  is empty, the program enters the main menue.
C
       open(97,FILE='scatbat')
       rewind(97)
       read(97,1000,ERR=1,END=1) la
C
C         Batchfile exists, call batchsubroutine
C
          if (la.eq.'begin') then
             batch=.true.
             call intlin
             goto 9999
          endif
C
C
C----------------------------------------------------------
C  Main menue section.
C
 1       call dispmm(n)
            if      (n.eq.1)   then
C               Generates Configurations
                call gencon
            elseif  (n.eq.2)   then
C               Calculates Bornamplitude
                call sborn
            elseif (n.eq.3)    then
C               s-wave multiple scattering
                call smult
            elseif (n.eq.4)    then
C               Enter command line interpreter
                call intlin
            elseif (n.eq.7)    then
C               Change Options+Files
                call option
            elseif  (n.eq.8)   then
C               Edit Parameters
                call edpar
            elseif (n.eq.9)    then
C               Save parameters
                call saveop
            elseif (n.eq.0)    then
C               Exit
                goto 9999
            else
                print *,'Illegal'
            endif
         goto 1
C
C----------------------------------------------------------
C  Close result file
C
 1000  format(a80)
 9999  close(99)
       end


       subroutine intlin
C
C  This subroutine is a primitive command line
C  interpreter. It can be used as alternative
C  to the MENUE environment. It is also used
C  for batch operation.
C  Input is taken from getlin. If batchmode is on
C  this means, that file 97 is read.
C
C  Auxilliary variables:
C
       real r
       integer i
       complex c
C
C  Input string:
C
       character*(128) ch
C
C  Externals
C
       integer rdint
       real    rdreal
       complex rdcmpl
C
C  Globals
C
C  This is the maximum number of particles
       parameter(maxpnm=20)
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C----------------------------------------------------------
C  Central loop
C
 100   continue
C
C Readline
C
         if (.not.batch) then
            write(*,*) ':'
         endif
         call getlin(ch)
         if (ch.eq.'#END') then
C             End of inputfile reached
              goto 9999
         elseif (ch.eq.'#ERR') then
C             Read error
              call error
              stop
         elseif (ch.eq.'config') then
C             Make configurations
              call gencon
         elseif (ch.eq.'smult')  then
C             s-wave multiple scattering
              call smult
         elseif (ch.eq.'sborn')   then
C             calculate bornamplitude
              call sborn
         elseif (ch.eq.'pborn')   then
C             calculate bornamplitude
              call pborn
         elseif (ch.eq.'f2')   then
C             calculate  p wave 2nd order amplitude
              call f2
C
C   Parameters can be mofified here
C
         elseif (ch.eq.'setk')   then
              r=rdreal(ch)
              if (ch.eq.'#OK')   then
                 k=r
              else
                 call error
              endif
          elseif (ch.eq.'seth')   then
              r=rdreal(ch)
              if (ch.eq.'#OK')   then
                 hardc=r
              else
                 call error
              endif
          elseif (ch.eq.'setnc')   then
              i=rdint(ch)
              if (ch.eq.'#OK')   then
                 nconf=i
              else
                 call error
              endif
          elseif (ch.eq.'setpn')   then
              i=rdint(ch)
              if (ch.eq.'#OK')   then
                 pnum=i
              else
                 call error
              endif
          elseif (ch.eq.'setth')   then
              i=rdint(ch)
              if ((ch.eq.'#OK').and.(thetas.le.180)) then
                 thetas=i
              else
                 call error
              endif
         elseif (ch.eq.'setw')   then
              r=rdreal(ch)
              if (ch.eq.'#OK')   then
                 width=r
              else
                 call error
              endif
         elseif (ch.eq.'setf0')   then
              c=rdcmpl(ch)
              if (ch.eq.'#OK')   then
                 f0=c
              else
                 call error
              endif
         elseif (ch.eq.'setf1')   then
              c=rdcmpl(ch)
              if (ch.eq.'#OK')   then
                 f1=c
              else
                 call error
              endif
C
C  End of parameter modification section
C
         elseif (ch.eq.'menue') then
C             Terminate commandline interpreter
              goto 9999
         elseif (ch.eq.'end') then
C             End program
              close(99)
              stop
         endif
C      Next line
 200   continue
       goto 100

9999   end


       subroutine error
C
C  Error message subroutine of intlin
C
       logical batch,makprt
       common/flags/makprt,batch
C
       if (makprt) then
       endif
       if (batch) then
           write(99,*) 'INPUT ERROR IN BATCHFILE'
           close(99)
           stop
       else
           print *,'INPUT ERROR'
       endif
       end


       integer function rdint(ch)
C
C This function read an integer from
C the input stream  getlin
C The variable ch is used here as input variable
C AND
C as flag for the calling procedure
C only if ch='#OK' after return, the result is valid !!!
C
       character*(*)     ch
       integer i
C
1      format(i8)
C
C Read line
C
       call getlin(ch)
       if ((ch.eq.'#ERR').or.(ch.eq.'#END')) then
          goto 9999
       else
          read(ch,1,ERR=100) i
          rdint=i
          ch='#OK'
          goto 9999
       endif
C
C Integer error trap
C
100    continue
       ch='#INTERR'
       rdint=0
C
9999   end

       real function rdreal(ch)
C
C This function read a real number from
C the input stream  getlin
C The variable ch is used here as input variable
C AND
C as flag for the calling procedure
C only if ch='#OK' after return, the result is valid !!!
C
       character*(*)   ch
       real r
C
1      format(g14.8)
C
C Read line
C
       call getlin(ch)
       if ((ch.eq.'#ERR').or.(ch.eq.'#END')) then
          rdreal=0.0
          goto 9999
       else
          read(ch,1,ERR=100) r
          rdreal=r
          ch='#OK'
          goto 9999
       endif
C
C Real error trap
C
100    continue
       ch='#REALERR'
       rdreal=0.0
C
9999   end




       complex function rdcmpl(ch)
C
C This function read a complex number from
C the input stream  getlin
C The variable ch is used here as input variable
C AND
C as flag for the calling procedure
C only if ch='#OK' after return, the result is valid !!!
C
       character*(*) ch
       real r1,r2
       complex c
C
1     format(2g14.8)
C
C Read line
C
       call getlin(ch)
       if ((ch.eq.'#ERR').or.(ch.eq.'#END')) then
          rdcmpl=(0.0,0.0)
          goto 9999
       else
          read(ch,1,ERR=100) r1,r2
          rdcmpl=r1+(0.0,1.0)*r2
          ch='#OK'
          goto 9999
       endif
C
C Complex error trap
C
100    continue
       ch='#CMPLERR'
       rdcmpl=(0.0,0.0)
C
9999   end


       subroutine getlin(ch)
C
C  This subroutine reads a line from an input divice
C  if batch is on, the input divice is the file 97
C  else the line is read from the console
C  In batchmode the line ch is always written to the
C  resultfile
C
       character*(*) ch
C
C  Globals
C
       logical makprt,batch
       common/flags/makprt,batch

       if (batch) then
           read(97,100,END=110,ERR=120) ch
           write(99,*) ch
       else
           read(*,100,END=110,ERR=120) ch
       endif
       goto 9999
C
C END OF FILE
C
 110   ch='#END'
       goto 9999
C
C READ ERROR
C
 120   ch='#ERR'
       goto 9999

 100   format(a80)
9999   end



       subroutine dispmm(n)
C
C      This subroutine displays the main menue items
C      and returns the result in n
C
       integer n
C
       print *,'   MAIN MENUE:                       '
       print *,'                                     '
       print *,'   (1)     Generate Configurations   '
       print *,'   (2)     Calculate Bornamplitude   '
       print *,'   (3)     s-wave Mult. Scatt.       '
       print *,'   (4)     Enter DM                  '
       print *,'   (7)     Options+Files             '
       print *,'   (8)     Edit Parameters           '
       print *,'   (9)     Save parameters           '
       print *,'   (0)     End Program               '
       print *,'                                     '
       print *,'?'
       read *,n
C
       end


       subroutine gencon
C
C   This subroutine generates configurations for
C   stochastical integration. The configurations are
C   gaussian ( EXP(-(X*X)/(width*width) ).
C   A hardcoreradius can be included (option!)
C   The Configurations are stored in the file config.
C   The format of the file is
C       titel  (this program generates the title 'gaussian')
C       pnum   (Number of scatterers)
C       nconf  (Number of configurations)
C       width
C       hardc
C       x1   y1   z1
C       . . .
C

C  This is the maximum number of particles
       parameter(maxpnm=20)
C
C   Auxilliary variables of the subroutine
C
C   Counters etc
       integer i,j,l,j1,j2
       real d,y,h2
C   maxpnm 3-d vectors
       real x(3,maxpnm)
C
C
C Commonblock from Main
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,th
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C      Initialisation
C
       h2=hardc*hardc
C
C      Warning if file exists
C
C not implemented
C
C      Open file and print initial values
C
       open(1,FILE=config)
       rewind(1)
       write(1,*) 'gaussian'
       write(1,100) pnum
       write(1,100) nconf
       write(1,*) width
       write(1,*) hardc
100    format(i8)
C
C  Print initial message
C
       if (.NOT.batch) then
          print *,' '
          print *,nconf,' gaussian configurations of '
          print *,pnum,' particles will be generated.'
          print *,'Width           : ',width
          print *,'Hardcore radius : ',hardc
          print *,' '
       endif
C
C      Generate configurations
C
C      nconf configurations
       do 1 i=1,nconf
C
C         of pnum particles
C
1000      continue
          do 2 j=1,pnum
C
C            in 3-d
C
             do 2 l=1,3
                x(l,j)=gauss(width)
 3           continue
 2        continue
C
C         Hardcoretest
C
          if (hardc.ne.0.0) then
             do 1001 j1=1,pnum
                do 1001 j2=(j1+1),pnum
C Distance between scatterer j1 and j2
                   d=0
                   do 1002 l=1,3
                      y=x(l,j1)-x(l,j2)
                      d=d+y*y
1002               continue
C Test
                   if (h2.gt.d) then
C                     Reject configuration
                      goto 1000
                   endif
1001         continue
          endif
C
C         Print configuration
C
          do 4 j=1,pnum
             write(1,101) x(1,j),x(2,j),x(3,j)
 101      format(3e12.4)
 4        continue
C
C      next configuration
C
 1     continue
       close(1)
C
C   Print message
C
       if (.NOT.batch)  then
          print *,'Information is stored in file ',config
          print *,' '
       endif
C
C   Print message to result file
C   if batchmode is on
C
       if (batch) then
         write(99,*) '--------------------------------------------'
         write(99,*) 'Subroutine gencon :'
         write(99,*) nconf,' gaussian configurations of '
         write(99,*) pnum,' particles generated.'
         write(99,*) 'Width           : ',width
         write(99,*) 'Hardcore radius : ',hardc
         write(99,*) 'Information is stored in file ',config
         write(99,*) '--------------------------------------------'
       endif

       end

       real function gauss(a)
C
C   Produces a gaussian
C   EXP(-x*x/(width*width))
C   random number
C   (very ineffective!)
C
       gauss=a*sqrt(-log(ran()))*cos(2*3.1415*ran())

       end

       REAL FUNCTION RAN()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE CALCULATES A RANDOM NUMBER                             C
C                                                                      C
C  It is taken from Dieter Stolls programm MUSH                        C                                            C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       DOUBLE PRECISION X,XK,PM,A,C
       PARAMETER(A=14662125.D00,C=13136923.D00,PM=2.D00**48)

       COMMON/RG/X,XK

       XK=XK+.5D00
       X=A*X+DINT(C*XK)
       X=X-DINT(X/PM)*PM
       RAN=REAL(X/PM)

       END


       subroutine sborn
C
C  This routine calculates the Born amplitude of
C  a given distribution of scatters.
C  The configurations are take from config.
C  The angular steps will be taken according
C  to the predefined values in the main program
C  If the distributian is gaussian, the exact result
C  is also written to the result file.

C
C  Auxilliary variables
C
C  Counters etc.
       integer i,j,l,t
       real qa
       character ch
C  Title of config
       character*(20) title
C  A vector
       real x(3)
C  Summation variable
       complex sum(0:180)
C
C  Momentum vectors
C  Maximum 180 vectors = 1 degree step
C
       real ki(3)
       real ko(3,0:180)
       real q(3,0:180)
C
C  Bornamplitude
C  fb  : stochastic result
C  fbt : analytic result
C
       complex  fb(0:180),fbt(0:180)
C
C  The following variables exist also as globals.
C  These are there local counterparts.
C  The loaded config file cannot change the global parameters.
C
       integer nconfl
       integer pnuml
       real widthl,hardcl
C
C  Globals as usual
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C  Start message
C
       if (.NOT.batch) then
          print *,' '
          print *,' Beginning s-wave Born calculation.'
          print *,' Momentum     : ',k
          print *,' Scat. Ampl.  : ',f0
       endif
C
C  Generate in and out  momentum vectors
C
C  The incoming momentum is in z direction
C
       ki(1)=0.0
       ki(2)=0.0
       ki(3)=k
C
C  The outgoing momentum is in the y-z plane
C
       do 10 i=0,thetas
          ko(1,i)=0.
          ko(2,i)=k*sin(pi/thetas*i)
          ko(3,i)=k*cos(pi/thetas*i)
 10    continue
C
C  Here ki-ko=q is calculated.
C  Mind the sign, q is incoming minus outgoing !
C
       do 20 i=0,thetas
          q(1,i)=0.0
          q(2,i)=ko(2,i)
          q(3,i)=k-ko(3,i)
 20    continue
C
C  Open config file and load parameters from it.
C
       open(1,FILE=config)
       rewind(1)
       read(1,'(a20)') title
       read(1,*)     pnuml
       read(1,*)     nconfl
       read(1,*)     widthl
       read(1,*)     hardcl
C
C  Print message  ( in interactive mode )
C
       if (.NOT.batch) then
          print *,' Opened configuration file :',config
          print *,' Distribution : ',title
          print *,' Particles    : ',pnuml
          print *,' Config.      : ',nconfl
          print *,' Width        : ',widthl
          print *,' Hardcore     : ',hardcl
        endif
C
C  Begin calculation
C
C Set summation to zero
       do 21 t=0,20
         sum(t)=(0.0,0.0)
 21    continue
C
C      Configurations
C
       rms=0.0
       do 30 i=1,nconfl
C         Particles
          do 30 j=1,pnuml
C            Read vector
             read(1,*) x(1),x(2),x(3)
               rms=rms+x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
C              Scattering angles
               do 25 t=0,thetas
                 sum(t)=sum(t)+
     *           exp((0.0,1.0)*(q(1,t)*x(1)+q(2,t)*x(2)+q(3,t)*x(3)))
 25           continue
C
C   Next particle, next configuration
C
 30    continue
C
C   Store result
C
       rms=rms/nconfl/pnuml
       print *,rms
       write(99,*)'RMS=',rms

       do  35 t=0,thetas
          fb(t)=f0*sum(t)/nconfl
C
C   Calculate analytic result
C
C   Calculate qa from k and theta
          qa=k*k*2*(1-cos(pi*t/thetas))
C   Scattering amplitude
          fbt(t)=f0*pnuml*exp(-qa*widthl*widthl/4.0)
C
C  Next angle
C
 35    continue
C
C   Output results
C
C   Resultfile:
          write(99,*)' ----------------------------------------'
          write(99,*)' Beginning Born calculation.'
          write(99,*)' Momentum    : ',k
          write(99,*)' Scat. ampl. : ',f0
          write(99,*)' Opened configuration file : ',config
          write(99,*)' Distribution : ',title
          write(99,*)' Particles    : ',pnuml
          write(99,*)' Config.      : ',nconfl
          write(99,*)' Width        : ',widthl
          write(99,*)' Hardcore     : ',hardcl
          write(99,*)
     *     ' Theta               Stochastic              Analytic '
          do 40 t=0,thetas
             write(99,*) (180/thetas*t),fb(t),fbt(t)
 40       continue
          write(99,*)'----------------------------------------'
C  Console
       if (.NOT.batch) then
          print *,'Theta         Stochastic          Analytic'
          do 50 t=0,thetas
             print *, (180/thetas*t),fb(t),fbt(t)
 50       continue
       print *,'Enter c to continue'
       read(*,1000,ERR=9999) ch
1000   format(a1)
       endif

9999   end


       subroutine edpar
C
C Here the user can change the global parameters
C (some of them).
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C Auxilliary variables of the subroutine
C
       integer n
C
C Display values
C
 1     print *,' '
       print *,'  PARAMETERS :'
       print *,' '
       print *,'(1)   Number of particles :  ',pnum
       print *,'(2)   Momentum            :  ',k
       print *,'(3)   Scatt. Ampl. f0     :  ',f0
       print *,'(4)   Scatt. Ampl. f1     :  ',f1
       print *,'(5)   Width               :  ',width
       print *,'(6)   Hardcore radius     :  ',hardc
       print *,'(7)   Numb. of Conf.      :  ',nconf
       print *,' '
       print *,'(0)   Main menue              '
       print *,' '
       print *,' ?'
       read *,n

       if (n.eq.1) then
          print *,'Number of particles :  '
          read  *,pnum
       elseif  (n.eq.2) then
          print *,'Momentum :  '
          read  *,k
       elseif  (n.eq.3) then
          print *,'s-wave scatt. Ampl. :  '
          read *,f0
       elseif  (n.eq.4) then
          print *,'p-wave scatt. Ampl. :  '
          read *,f1
       elseif  (n.eq.5) then
          print *,'Width :  '
          read *,width
       elseif  (n.eq.6) then
          print *,'Hardcore radius :  '
          read *,hardc
       elseif  (n.eq.7) then
          print *,'Number of configurations :  '
          read *,nconf
      elseif (n.eq.0) then
          goto 99
       else
          print *,'Illegal'
       endif
       goto 1

99     end


       subroutine option
C
C Here filenames and flag are controlled.
C The structure is the same as in edpar
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/flags/makprt,batch
C
C Auxilliary variables of the subroutine
C
       integer n
C
C      External functions
C
       character*(4) onoff
C
C Display values
C
 1     print *,' '
       print *,'  PARAMETERS :'
       print *,' '
       print *,'(1)   Configurationfile   :  ',config
       print *,'(2)   Pfile               :  ',prtfil
       print *,'(3)   P                   :  ',onoff(makprt)
       print *,' '
       print *,'(0)   Main menue              '
       print *,' '
       print *,' ?'
       read *,n


       if  (n.eq.1) then
          print *,'Configurationfile  :  '
          read  10,config
       elseif  (n.eq.2) then
          print *,'Pfile  :  '
          read  10,prtfil
       elseif  (n.eq.3) then
          makprt=(.NOT.makprt)
       elseif (n.eq.0) then
          goto 99
       else
          print *,'Illegal'
       endif
       goto 1


 10    format(a20)
 99    end


       character*(4) function onoff(l)
       logical l
       if (l) then
          onoff='ON'
       else
          onoff='OFF'
       endif
       end


        subroutine saveop
C
C    Here, the file scatpar is generated.
C    It contains a set of parameters, which is loaded
C    at program start
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch


          open(98,FILE='scatpar')
          rewind(98)
          write(98,*) pnum
          write(98,202) config
          write(98,202) prtfil
          write(98,*) makprt
          write(98,*) k
          write(98,*) f0
          write(98,*) f1
          write(98,*) width
          write(98,*) hardc
          write(98,*) nconf
 202      format(a20)
          close(98)
       end


       subroutine smult
C
C Here the s-wave multiple scattering amplitude is calculated
C
C Maximum number of particles
       parameter(maxpnm=20)
C
C  Auxilliary variables
C
C  Counters etc.
       integer i,j,l,t,ai,aj,bi,bj
       real qa,d,y
       character ch
C  Title of config
       character*(20) title
C  A configuration
       real x(3,maxpnm)
C  Summation variable
       complex sum
C  Greensfunction matrix
       complex a(maxpnm,maxpnm)
C  State vectors
       complex b(maxpnm),ci(maxpnm),co(maxpnm)
C Auxilliary variable for Gauss algorithm
       complex h(maxpnm),hs
C
C  Momentum vectors
C  Maximum 180 vectors = 1 degree step
C
       real ki(3)
       real ko(3,0:180)
       real q(3,0:180)
C
C
C  fs  : stochastic result
C
       complex  fs(0:180)
C
C  The following variables exist also as globals.
C  These are there local counterparts.
C  The loaded config file cannot change the global parameters.
C
       integer nconfl
       integer pnuml
       real widthl,hardcl
C
C  Globals as usual
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C  Start message
C
       if (.NOT.batch) then
          print *,' '
          print *,' Beginning s-wave multiple scattering calculation.'
          print *,' Momentum     : ',k
          print *,' Scat. Ampl.  : ',f0
       endif
C
C  Generate in and out  momentum vectors
C
C  The incoming momentum is in z direction
C
       ki(1)=0.0
       ki(2)=0.0
       ki(3)=k
C
C  The outgoing momentum is in the y-z plane
C
       do 10 i=0,thetas
          ko(1,i)=0.
          ko(2,i)=k*sin(pi/thetas*i)
          ko(3,i)=k*cos(pi/thetas*i)
 10    continue
C
C  Here ki-ko=q is calculated.
C  Mind the sign, q is incoming minus outgoing !
C
       do 20 i=0,thetas
          q(1,i)=0.0
          q(2,i)=ko(2,i)
          q(3,i)=k-ko(3,i)
 20    continue
C
C  Open config file and load parameters from it.
C
       open(1,FILE=config)
       rewind(1)
       read(1,'(a20)') title
       read(1,*)     pnuml
       read(1,*)     nconfl
       read(1,*)     widthl
       read(1,*)     hardcl
C
C  Print message  ( in interactive mode )
C
       if (.NOT.batch) then
          print *,' Opened configuration file :',config
          print *,' Distribution : ',title
          print *,' Particles    : ',pnuml
          print *,' Config.      : ',nconfl
          print *,' Width        : ',widthl
          print *,' Hardcore     : ',hardcl
        endif
C
C  Begin calculation
C
C Set fs to zero
       do 21 t=0,20
         fs(t)=(0.0,0.0)
 21    continue
C
C      Configurations
C
       do 30 i=1,nconfl
C         Particles
          do 31 j=1,pnuml
C            Read vectors
             read(1,*) x(1,j),x(2,j),x(3,j)
 31       continue
C
C  Solve multiple scattering equations for the given
C  Configuration  x(i,j)
C
C First: Build up Greensfunction matrix
C
          do 100 ai=1,pnuml
            do 100 aj=ai,pnuml
              if (ai.eq.aj) then
                 a(ai,aj)=1.0
              else
C Distance between scatterer ai and aj
                 d=0
                 do 110 l=1,3
                    y=x(l,ai)-x(l,aj)
                    d=d+y*y
 110             continue
                 d=sqrt(d)
C Greensfunction  * f0
                 a(ai,aj)=f0*exp((0.0,1.0)*k*d)/d
                 a(aj,ai)=a(ai,aj)
              endif
 100      continue

C
C  Second : Construct incoming state vector
C
          do 120 bi=1,pnuml
             ci(bi)=exp((0.0,1.0)*
     *             ( ki(1)*x(1,bi)+
     *               ki(2)*x(2,bi)+
     *               ki(3)*x(3,bi) ))
 120      continue
C
C  Third : Solve linear equations
C  Solution : b
C  Method : Bronstein, Taschenbuch der Mathematik
C           S.735 ( Gaussches Eliminationsverfahren )
C
C Transform A to triangular shape
          do 200 ai=1,(pnuml-1)
             do 230 aj=(ai+1),pnuml
                h(aj)=-a(aj,ai)/a(ai,ai)
                do 220 bi=1,pnuml
                  if (bi.eq.ai) then
                     a(aj,bi)=0.0
                  else
                     a(aj,bi)=a(aj,bi)+h(aj)*a(ai,bi)
                  endif
 220            continue
                ci(aj)=ci(aj)+h(aj)*ci(ai)
 230         continue
 200      continue
C Generate solution vector
          b(pnuml)=ci(pnuml)/a(pnuml,pnuml)
          do 300 ai=(pnuml-1),1,-1
             hs=0.0
             do 310 aj=1,(pnuml-ai)
                hs=hs+a(ai,ai+aj)*b(ai+aj)
 310         continue
          b(ai)=(ci(ai)-hs)/a(ai,ai)
 300      continue

C
C  Fourth : Project solution on outgoing vector
C  Note this is the first step, where the angular
C  dependence comes in
C
C  Scattering angles:
          do 130 t=0,thetas
C
C  Outgoing state vector for angle t
C
            do 140 bi=1,pnuml
              co(bi)=exp((0.0,-1.0)*
     *             ( ko(1,t)*x(1,bi)+
     *               ko(2,t)*x(2,bi)+
     *               ko(3,t)*x(3,bi) ))
 140        continue
C
C  Projection
C
            sum=0
            do 150 bi=1,pnuml
                sum=sum+co(bi)*b(bi)
 150        continue
C
C   Add individual configuration
C
            fs(t)=fs(t)+(-1.0)*f0*sum
C
C   Next angle
C
 130      continue
C
C   Next configuration
C
 30    continue
C
C    Divide result
C
       do 400 t=0,thetas
         fs(t)=fs(t)/nconfl
 400   continue
C
C   Output results
C
C   Resultfile:
          write(99,*)' ----------------------------------------'
          write(99,*)' Beginning s-wave ms calculation.'
          write(99,*)' Momentum    : ',k
          write(99,*)' Scat. ampl. : ',f0
          write(99,*)' Opened configuration file : ',config
          write(99,*)' Distribution : ',title
          write(99,*)' Particles    : ',pnuml
          write(99,*)' Config.      : ',nconfl
          write(99,*)' Width        : ',widthl
          write(99,*)' Hardcore     : ',hardcl
          write(99,*)
     *     ' Theta               fs '
          do 40 t=0,thetas
             write(99,*) (180/thetas*t),fs(t)
 40       continue
          write(99,*)'----------------------------------------'
C  Console
       if (.NOT.batch) then
          print *,'Theta           fs'
          do 50 t=0,thetas
             print *, (180/thetas*t),fs(t)
 50       continue
         print *,'Enter c to continue'
         read(*,*,ERR=9999) ch
       endif
9999   end



       subroutine f2
C
C  This routine calculates the pure p-wave
C  second order born tern amplitude
C  (in forward direction).
C  Motivation : Thesis Helgert inequality of integeral
C  I1-I3
C
C  Extension: The rescattering 3rd order term
C             is included to check the problems
C             of P200 ff
C
C  The configurations are take from config.
C
       parameter(maxpnm=20)

C
C  Auxilliary variables
C
C  Counters etc.
       integer i,j,l,t,n,m
       real d,y,dz,kd,yd,sq,sq3,dx,dy,ay2,ay0
C  Title of config
       character*(20) title
C  Vectors
       real x(3,maxpnm)
C  Summation variable
       complex sum1,sum2,funreg,freg
       complex sum3,sum4,kiss,fa(10),h0,h2,sum5,sum6
C
C
C  These are there local counterparts of the global variables.
C  The loaded config file cannot change the global parameters.
C
       integer nconfl
       integer pnuml
       real widthl,hardcl
C
C  Globals as usual
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C  Initialisation
C
       sq=sqrt(2.0)
       sq3=sq*sq*sq

       fa(1)=(.0010,.0010)
       fa(2)=(.002,.002)
       fa(3)=(.0030,.0030)
       fa(4)=(.004,.004)
       fa(5)=(.005,.0050)
       fa(6)=(.006,.006)
       fa(7)=(.0070,.0070)
       fa(8)=(.008,.008)
       fa(9)=(.0090,.0090)
       fa(10)=(.0100,.0100)
C
C  Start message
C
       if (.NOT.batch) then
          print *,' '
          print *,' Beginning p-wave 2nd order calculation.'
          print *,' Momentum     : ',k
          print *,' Scat. Ampl.  : ',f1
       endif
C
C  Open config file and load parameters from it.
C
       open(1,FILE=config)
       rewind(1)
       read(1,'(a20)') title
       read(1,*)     pnuml
       read(1,*)     nconfl
       read(1,*)     widthl
       read(1,*)     hardcl
C
C  Print message  ( in interactive mode )
C
       if (.NOT.batch) then
          print *,' Opened configuration file :',config
          print *,' Distribution : ',title
          print *,' Particles    : ',pnuml
          print *,' Config.      : ',nconfl
          print *,' Width        : ',widthl
          print *,' Hardcore     : ',hardcl
        endif
C
C  Begin calculation
C
C Set summation to zero
         sum1=0.
         sum2=0.
         sum3=0.
         sum4=0.
         sum5=0.
         sum6=0.
C
C      Configurations
C
       do 40 i=1,nconfl
C         Particles
          do 30 j=1,pnuml
C            Read vectors
             read(1,*) x(1,j),x(2,j),x(3,j)
30        continue
C
C  All possible combinations of particles
C
          do 31 n=1,pnuml
            do 31 m=1,pnuml
            if (m.ne.n) then
C Distance between scatterer n and m
          d=0.
          do 110 l=1,3
             y=x(l,n)-x(l,m)
             d=d+y*y
 110      continue
          d=sqrt(d)
C cos(theta)
          dz=(x(3,n)-x(3,m))/d
          dy=(x(2,n)-x(2,m))/d
          dx=(x(1,n)-x(1,m))/d
          kd=d*k
C
C  hankels and sper. harmonics
C
          h0=exp((0.0,1.0)*kd)/kd
          h2=h0*(3/(kd*kd)-1-(0.0,3.0)/kd)
          ay0=(dx*dx+dy*dy)*dz*dz
          ay2=(2.*dz*dz-dx*dx-dy*dy)
C
C Unregulated scattering amplitude
C
          funreg=
     A      exp((0.0,1.0)*kd*(-dz+1))/kd*
     B      (
     C        1.0-
     D        (3.0/(kd*kd)-(0.0,3.0)/kd-1.0)*(3.0*dz*dz-1.0)
     E      )
C
C    third order backscattering
C
          freg=h2*h2*(9*ay0+ay2*ay2)+
     A         h0*h2*(-2*ay2)+
     B         h0*h0
C
C   next particle
C
          sum1=sum1+funreg
          sum5=sum5+freg
           endif
31         continue


C
C   Next configuration
C
40      continue
        do 50 i=1,10
          f1=fa(i)
C
C   Store result
C
       sum2=3*f1*f1*k*sum1/nconfl+3*f1*pnuml
       sum6=3*f1*f1*f1*k*k*sum5/nconfl+sum2
       sum4=sum2+3*sqrt(2./pi)/(widthl*widthl*widthl)*f1*f1/(k*k)
     A              *pnuml*(pnuml-1)
C
C   Output results
C
C   Resultfile:
C
        write(99,*)' ----------------------------------------'
        write(99,*)' p-wave 2nd order calculation.'
        write(99,*)' Momentum    : ',k
        write(99,*)' Scat. ampl. : ',f1
        write(99,*)' Opened configuration file : ',config
        write(99,*)' Distribution : ',title
        write(99,*)' Particles    : ',pnuml
        write(99,*)' Config.      : ',nconfl
        write(99,*)' Width        : ',widthl
        write(99,*)' Hardcore     : ',hardcl
        write(99,*)
        write(99,*)
        write(99,*)' Scattering Amplitude:'
        write(99,*)'Amplitude   :  ',sum2
        write(99,*)'Zusatz      :  ',sum4
        write(99,*)'Rueckstreu  :  ',sum6

        write(96,8888) cplo,cabs(f1),aimag(sum2)
        write(96,8888) cplo+1,cabs(f1),real(sum2)
        write(96,8888) cplo+2,cabs(f1),aimag(sum4)
        write(96,8888) cplo+3,cabs(f1),real(sum4)
        write(96,8888) cplo+4,cabs(f1),aimag(sum6)
        write(96,8888) cplo+5,cabs(f1),real(sum6)


8888    format(i3,2e12.4)
C
C   Standard output
C
        print *,' Scattering Amplitude:'
        print *,'Amplitude   :  ',sum2
        print *,'+Zusatz     :  ',sum4
        print *,'+Ruechstreu :  ',sum6
50     continue

       cplo=cplo+6
9999   end



       subroutine pborn
C
C  This routine calculates the Born amplitude of
C  a given distribution of scatters for
C  PURE P-WAVE SCATTERING up to second order.
C  The configurations are take from config.
C  The angular steps will be taken according
C  to the predefined values in the main program
C  If the distributian is gaussian, the
C  first order exact result
C  is also written to the result file.
C
       parameter(maxpnm=20)

C
C  Auxilliary variables
C
C  Counters etc.
       integer i,j,l,t,n,m
       real qa,da,kd
       character*(4) ch
       complex e1,e2
C  Title of config
       character*(20) title
C  A vector
       real x(3,maxpnm)
       real d(3)
C  Summation variable
       complex sum1(0:180),sum2(0:180)
C
C  Momentum vectors
C  Maximum 180 vectors = 1 degree step
C
       real ki(3)
       real ko(3,0:180)
       real q(3,0:180)
C
C  Bornamplitude
C  fb1  : stochastic result first order
C  fb2  : stochastic result second order
C  fbt1 : analytic result
C
       complex  fb1(0:180),fbt1(0:180),fb2(0:180)
C
C  The following variables exist also as globals.
C  These are there local counterparts.
C  The loaded config file cannot change the global parameters.
C
       integer nconfl
       integer pnuml
       real widthl,hardcl
C
C  Globals as usual
C
       integer   pnum
       character*(20) config
       logical   makprt,batch
       character*(20) prtfil
       real      k
       real      width,hardc
       complex f0,f1
       integer   nconf
       character*(20) result
       character*(20) plotf
       integer cplo
       integer thetas
       real pi

       common/filnam/config,prtfil,result,plotf
       common/ipar/pnum,nconf,cplo,thetas
       common/rpar/k,width,hardc,pi
       common/cpar/f0,f1
       common/flags/makprt,batch
C
C  Start message
C
       if (.NOT.batch) then
          print *,' '
          print *,' Beginning p-wave Born calculation.'
          print *,' Momentum     : ',k
          print *,' Scat. Ampl.  : ',f1
       endif
C
C  Generate in and out  momentum vectors
C
C  The incoming momentum is in z direction
C
       ki(1)=0.0
       ki(2)=0.0
       ki(3)=k
C
C  The outgoing momentum is in the y-z plane
C
       do 10 i=0,thetas
          ko(1,i)=0.
          ko(2,i)=k*sin(pi/thetas*i)
          ko(3,i)=k*cos(pi/thetas*i)
 10    continue
C
C  Here ki-ko=q is calculated.
C  Mind the sign, q is incoming minus outgoing !
C
       do 20 i=0,thetas
          q(1,i)=0.0
          q(2,i)=ko(2,i)
          q(3,i)=k-ko(3,i)
 20    continue
C
C  Open config file and load parameters from it.
C
       open(1,FILE=config)
       rewind(1)
       read(1,'(a20)') title
       read(1,*)     pnuml
       read(1,*)     nconfl
       read(1,*)     widthl
       read(1,*)     hardcl
C
C  Print message  ( in interactive mode )
C
       if (.NOT.batch) then
          print *,' Opened configuration file :',config
          print *,' Distribution : ',title
          print *,' Particles    : ',pnuml
          print *,' Config.      : ',nconfl
          print *,' Width        : ',widthl
          print *,' Hardcore     : ',hardcl
        endif
C
C  Begin calculation
C
C Set summation to zero
       do 21 t=0,thetas
         sum1(t)=(0.0,0.0)
         sum2(t)=(0.0,0.0)
 21    continue
C
C      Configurations
C
       do 30 i=1,nconfl
C         Particles
          do 31 j=1,pnuml
C            Read vector
             read(1,*) x(1,j),x(2,j),x(3,j)
31        continue
C
C  Scattering angles
C
C
C   First order term:
        do 25 j=1,pnuml
          do 25 t=0,thetas
           sum1(t)=sum1(t)+
     A     exp( (0.0,1.0)
     B           *(q(1,t)*x(1,j)+q(2,t)*x(2,j)+q(3,t)*x(3,j))
     C        )
 25       continue
C
C  Second order term:
C
        do 26 t=0,thetas
         do 26 n=1,pnuml
           do 26 m=1,pnuml
           if (n.ne.m) then
C
C   incomming and outgoing wave
C
            e1=
     A      exp( (0.0,-1.0)
     B            *(ko(1,t)*x(1,n)+ko(2,t)*x(2,n)+ko(3,t)*x(3,n))
     C          )
            e2=
     A      exp( (0.0,1.0)
     B            *(ki(1)*x(1,m)+ki(2)*x(2,m)+ki(3)*x(3,m))
     C          )
C
C   calculated distance
C
            da=0.0
            do 28 j=1,3
               d(j)=x(j,n)-x(j,m)
               da=da+d(j)*d(j)
 28         continue
            da=sqrt(da)
            kd=k*da
C
C   normalize d
C
            do 29 j=1,3
               d(j)=d(j)/da
 29         continue
C
C   scattering amplitude  (unregulated)
C
            sum2(t)=sum2(t)+
     A         e1*e2*
     C         (
     B           (-sqrt(3.0)*ko(2,t)/k*d(2)*d(3)
     C            -          ko(3,t)/k*(3*d(3)*d(3)-1) )
     D           *(3.0/(kd*kd)-(0.0,3.0)/kd-1)*exp((0.0,1.0)*kd)/kd
     E           +(ko(3,t)/k*exp((0.0,1.0)*kd)/kd)
     F         )
        endif
 26     continue

C
C   next configuration
C
 30    continue
C
C   Store result
C
       do  35 t=0,thetas
C
C   Calculate result
C
C   Calculate qa from k and theta
          qa=k*k*2*(1-cos(pi*t/thetas))
C   Scattering amplitude
          fb1(t)=3.0*f1*cos(pi*t/thetas)*sum1(t)/nconfl
          fbt1(t)=3*pnuml*f1*cos(pi*t/thetas)*
     A              exp(-qa/(widthl*widthl*2.0))
C  Factor 2 because only half of the dij are taken
          fb2(t)=3*f1*f1*k*sum2(t)/nconfl
C
C  Next angle
C
 35    continue
C
C   Output results
C
C   Resultfile:
          write(99,*)' ----------------------------------------'
          write(99,*)' Beginning p-Born calculation.'
          write(99,*)' Momentum    : ',k
          write(99,*)' Scat. ampl. : ',f1
          write(99,*)' Opened configuration file : ',config
          write(99,*)' Distribution : ',title
          write(99,*)' Particles    : ',pnuml
          write(99,*)' Config.      : ',nconfl
          write(99,*)' Width        : ',widthl
          write(99,*)' Hardcore     : ',hardcl
          write(99,*)' First order result:'
C Resultfile
          write(99,*)
     *     ' Theta               Stochastic              Analytic '
          do 40 t=0,thetas
             write(99,*) (180/thetas*t),fb1(t),fbt1(t)
 40       continue
          write(99,*)
          write(99,*)'Pure second order result:'
          do 41 t=0,thetas
             write(99,*) (180/thetas*t),fb2(t)
 41       continue
          write(99,*)
          write(99,*)'Sum of first and second order:'
          write(99,*)'Plotfilenumbers:'
          write(99,*)'Real part  ',cplo
          write(99,*)'Im.  part  ',cplo+1
          write(99,*)'X-sect    ',cplo+2
          do 42 t=0,thetas
             write(99,*) (180/thetas*t),fb1(t)+fb2(t)
      write(96,7777) cplo,(180./thetas*t),real(fb1(t)+fb2(t))
      write(96,7777) cplo+1,(180./thetas*t),aimag(fb1(t)+fb2(t))
      write(96,7777) cplo+2,(180./thetas*t),cabs(fb1(t)+fb2(t))**2
 42       continue
          cplo=cplo+2
          write(99,*)'----------------------------------------'
C  Console
       if (.NOT.batch) then
          print *,'First order result:'
          print *,'Theta         Stochastic          Analytic'
          do 50 t=0,thetas
             print *, (180/thetas*t),fb1(t),fbt1(t)
 50       continue
          print *,'Enter c to continue'
          read(*,8888,ERR=9999) ch
          print *,' '
          print *,'Pure second order result:'
          do 51 t=0,thetas
             print *, (180/thetas*t),fb2(t)
 51       continue
          print *,'Enter c to continue'
          read(*,8888,ERR=9999) ch
          print *,' '
          print *,'Sum of first and second order:'

          do 52 t=0,thetas
             print *, (180/thetas*t),fb1(t)+fb2(t)
 52       continue
          print *,'Enter c to continue'
          read(*,8888,ERR=9999) ch
       endif

7777   format(i3,2e12.4)
8888   format(a2)
9999   end
