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
       close(10)
