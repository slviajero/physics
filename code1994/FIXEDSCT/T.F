        k10=0
        write(17,1000)98,q,real(fkkp(1,k1,kf+1))/summe
        write(17,1000)99,q,imag(fkkp(1,k1,kf+1))/summe
        do 230 k2=1,kangle,mthet
          k10=k10+1
          write(15,1035)kf+1,k2-1.,10.*cabs(fkkp(k10,k1,kf+1)/
     1    summe)**2.,sqrt(qangle(1,k10)-qangle(2,k10)**2./kf),0.
230     continue

        WRITE(2,1030)Q
        IF (K1.EQ.1) WRITE(2,1050)
        K10=0
	gint = 0
        DO 250 K2=1,KANGLE,MTHET
c
c elastic scattering
c
          K10=K10+1
            WRITE(2,1040)K2-1,FKKP(k10,K1,Kf+1)/SUMme,10.*
     1            (CABS(FKKP(k10,K1,Kf+1))/SUMme)**2 
	  gint=gint+10.*cabs(fkkp(k10,k1,kf+1)/summe)**2
     1           *sin((k2-1)/180.*pi)*mthet/180.*pi
250       continue
	write(2,*) 'Elastic cross section   [mb]',gint*2*pi
	write(2,*) 'Total cross section     [mb]',
     1               40*pi/q*imag(fkkp(k1,1,kf+1)/summe)
	write(2,*) 'Inelastic cross section [mb]',
     1               40*pi/q*imag(fkkp(k1,1,kf+1)/summe)-gint*2*pi
        write(2,*) 'Elementary cross section [mb]',
     1               40*pi/q*imag(pleng0)

        write(18,1100)1,q,gint*2*pi
	write(18,1100)2,q,40*pi/q*imag(fkkp(k1,1,kf+1)/summe)
	write(18,1100)3,q,40*pi/q*imag(fkkp(k1,1,kf+1)/summe)-gint*2*pi
        write(18,1100)4,q,40*pi/q*imag(pleng0)

