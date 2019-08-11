        program r

        open(1,file='i')
        open(2,file='o')
        n=4
 100    read(1,*,end=9999,err=999) x,y
        write(2,1000) n,x,y
        goto 100

1000    format(i3,2e12.4)

999     if(mod(n,2).eq.0)then
           n=n+1
        else
           n=n+9
        endif
        goto 100

9999    end
