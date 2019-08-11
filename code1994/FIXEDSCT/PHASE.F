      program phase
c
c this program extracts ImF(0) and ReF(0)
c as a function of k out of the 
c pndbl output file forwar and writes
c    k and ImF/ReF
c on a output file
c
      open(1, file='forwar')

100   read(1,1000,end=9999) n,x,y
      if (n.eq.96) then
	 read(1,1000,end=9999) n1,x1,y1
	 write(*,*) x,y1/y
      endif
      goto 100
1000  format(i2,2e12.4)
9999  end
      
