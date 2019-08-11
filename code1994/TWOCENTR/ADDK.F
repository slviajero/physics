       program addk 
       implicit none
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
c
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
c addiere Increment zu k
c
       open(10,file='parin')
       write(10,*) mt,mp
       write(10,*) hb
       write(10,*) omega
       write(10,*) pin+.1
       write(10,*) tau
       write(10,*) nt,np
       close(10)

       end
