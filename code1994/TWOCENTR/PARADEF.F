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
