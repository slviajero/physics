       program ausw
       integer nlamx
       parameter(nlmax=20000)
       real q(0:nlmax)
       integer nl
       real dt
       integer i,k
       integer ngg
       real x0,x1,x2,x
       real sum
       real kk,s
       integer n
       real pi
c
c einlesen        
c
       pi=4.*atan(1.) 

       open(10,file='test')
       read (10,*) nl,dt
       do 110 i=0,nl
          read(10,*) q(i)
110    continue
       close(10)
      
c  
      open(11,file='res')
c
c gleichgewichtswert       
c
       print *,' Abschaetzung fuer GG-wert:'
       ngg=nl/3               
       n=nl-ngg+1
c
c Momente
c
       x1=0.
       x2=0.
       do 120 i=ngg,nl
          x1=x1+q(i)          
          x2=x2+q(i)*q(i)
120    continue           
       x1=x1/n
       x2=x2/n
       print *,'Momente :'
       print *,'<x>   = ', x1  
       print *,'<x^2> = ', x2 , 1/2.
       s=x2-x1**2

       print *,'Fouriertrafo :'
       do 130 kk=0.,6.,.2
          sum=0.
          do 140 i=ngg,nl
             sum=sum+cos(kk*q(i))
 140      continue
          sum=sum/n
c          print *,kk,sum,exp(-kk**2/4)
c          write(11,*) kk,log(sum),log(exp(-kk**2/4))
           print *,kk,sum,1/(kk**2+1.)
           write(11,*) kk,log(sum),log(1/(kk**2+1.))
 130   continue

       sum=0.
       print *,'Normierung :'
       do 150 i=ngg,nl
          x=q(i)
          sum=sum+exp(-abs(x))*exp(-(x-x1)**2/(2.*s))
150    continue
       sum=sum/n
       x0=1/sum*sqrt(2*pi*s)
       print *,' x0  =',x0
            
       close(11)
       end
