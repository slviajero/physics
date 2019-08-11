          program g

	  real a(1000,1000)

100       continue
	  do i=1,1000 
	    do j=1,1000
		 a(i,j)=i*j
            end do
          end do
          goto 100
	   
          end
