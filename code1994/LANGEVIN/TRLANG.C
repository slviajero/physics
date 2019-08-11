/*
    Berechnet die Streulaenge nach der Trace Methode
	(Gelman+Spruch). Dieses Programm hat experimentellen Charakter.
*/

/*
   Einige von PCCs Dingen
*/
#include <stdio.h>
#include <math.h>

/*
   Einige Dimensionen
   NBMAX maximale Zahl der beta Schritte
   NDIM  Raumdimensionen
   NDIST Groesse des Feldes fuer die Verteilung der
	 Pfadpunkte
   rmax  maximaler Radius fuer die Pfadpunkte
*/
#define NBMAX   800
#define NDIM      3
#define NDIST   100
#define rmax    30.
#define pi        M_PI
#define TRUE      1
#define FALSE     0

/*
  Random numbers
*/
      
#define a 14662125.
#define c 13136923.
      
static double xr,xk,tt,pm;

/*
   Einige Globals
   m   : Projektilmasse
   b   : Weite des Potentials;    br  : Referenzpotential;
   v0  : Staerke des Potentials   v0r : Referenzpotential;

   dt  : Langevin Zeitschritt
   ds=sqrt(dt*2);
   nc  : Zahl der Pfade
   ng  : Gleichgewichtslauf 

   db  : Beta Zeitschritt
   nb  : Zahl der Betaschritte
*/
float v0,m,b ;
float v0r,br;
float dt,ds;
int nc,ng;
float db;
int nb;
/*
   Langevin Felder
   w  : Zufallszahlen
   dx : speichert den Zwischenschritt
   x  : der Pfad
*/
static float x[NDIM+1][NBMAX+1];
static float dx[NDIM+1][NBMAX+1];
static float w[NDIM+1][NBMAX+1];
/*
   Das RDIST Feld, die Verteilung der Pfade wird darin 
   bestimmt.
*/
static int phir[NDIST];

main()
{ /* Variablen                     */ 

     int i,j;
     int n;
     float act1,actref1;
     float vv,vcorr;
     float erg,erg1;
     float r;
     char actwr;
     char rdtwr;

     extern float vpot();
     extern float vref();
     extern double ran();
 
     FILE *result, *action, *rdist, *fopen();
/*
     diverse Flags
        actwr : Zeitverlauf protokolieren
        rdtwr : RDIST messen 
*/
     actwr=TRUE;
     rdtwr=TRUE;
/*  
     INITIALIZATION OF THE RANDOM NUMBER GENERATOR 
*/
      pm=ldexp(1.,48);
      xr=pi*1.e11;
      xk=0.;
      for(i=0;i<20;i++) tt=ran();
/*
   diverse Files
*/
     result = fopen("tratio","w");
     action = fopen("tact","w");
     rdist  = fopen("rdist","w");

/* 
   Simulationsparameter lesen    
*/
  /* (1) physikalische Parameter : */
     printf(" Weite des Potentials   : \n ");
     scanf("%f",&b);
     printf(" Staerke des Potentials : \n ");
     scanf("%f",&v0);
     printf(" Masse des Projektile   : \n ");
     scanf("%f",&m);
  /* (2) physikalische Zeit */
     printf(" Zahl der beta Schritte : \n ");
     scanf("%d",&nb);
     printf(" db                     : \n ");
     scanf("%f",&db);
  /* (3) Langevin Zeit */
     printf(" Zahl der Langevin Schritte : \n ");
     scanf("%d",&nc);
     printf(" Gleichgewichtsschritte     : \n ");
     scanf("%d",&ng);
     printf(" dt                     : \n ");
     scanf("%f",&dt);
     ds=sqrt(dt*2.);
     printf("Langevin Simulation von Tr : \n");
     printf(" Pfad db = %f   nb = %d \n ",db,nb);
     printf(" Lang nc = %d   ng = %d   dt = %f \n",nc,ng,dt);
     printf(" Pot  b  = %f   v0 = %f   m  = %f \n",b,v0,m);
/* 
   Pfad initialisieren  
*/
     for(i=1;i<=nb;i++){
	for(j=1;j<=NDIM;j++){
	    x[j][i]=0.;
         };
     }
/*
   phir initialisieren
*/   
     for(i=0;i<NDIST;i++) phir[i]=0 ; 

/* 
     Referenzpotential bestimmen 
*/
     br  = b*sqrt(5.);
     v0r = v0*3*sqrt(pi/250.) ;
/* 
     Simulationsbegin 
*/
     erg=0.;
     for(n=1;n<=(nc+ng);n++){
/*  
       Potentialanteil der Wirkung berechnen

         hier wird etwas getrickst:
	 es wird sowohl der Korrekturfaktor fuer die Wirkung
	 als auch der Potentialteil des Driftterms berechnet.
	 In das Array dx wird die Potentialdrift geschrieben.
*/
	act1    = 0.;
	actref1 = 0.;
        for(i=1;i<=nb;i++){
/*        
            Radius und Potential ausrechnen 
*/
            r=0.;
            for(j=1;j<=NDIM;j++) r=r+x[j][i]*x[j][i];
            r=sqrt(r);
/*
   Verteilung der r bestimmen
*/
            if (rdtwr && i==1 && r<rmax) phir[(int)(r/rmax*NDIST)]++ ;
            vv=vpot(r)*db;
/* 
            (1) Potentialdrift 
*/
            for(j=1;j<=NDIM;j++) dx[j][i]=-vv/(b*b)*x[j][i]*dt;
            act1    = vv + act1    ;
            actref1 = vref(r) + actref1 ;
          }
          vcorr=-exp(-act1)/(exp(-act1)-1.);
          actref1 = actref1*db; /* Potentialteil der Wirkung fertig */
/*
           Ergebniss eintragen 
*/
          if (n>ng) {
                erg1 = (exp(-actref1)-1.)/(exp(-act1)-1.);
                erg  = erg+erg1;
/*
                falls actwr gesetzt: jeden 100. Schritt rausschreiben
*/
                if (actwr)
                    if ((n%100)==0) 
                       fprintf(action," %f %f \n",(n-ng)*dt,erg/(n-ng));
	      };
/* 
          (2) kinetische Energie + noise 
*/
          makran();
          for(i=2;i<nb;i++)
             for(j=1;j<=NDIM;j++)
	    	 dx[j][i]=m/db*(x[j][i+1]-2.*x[j][i]+x[j][i-1])*dt+
			          dx[j][i]*vcorr+w[j][i];
          for(j=1;j<=NDIM;j++)
             dx[j][1]  = m/db*(x[j][2]-2.*x[j][1]+x[j][nb])*dt+
			             dx[j][1]*vcorr+w[j][1];
          for(j=1;j<=NDIM;j++)
             dx[j][nb] = m/db*(x[j][1]-2.*x[j][nb]+x[j][nb-1])*dt+
	  		             dx[j][nb]*vcorr+w[j][nb];
/*
          Update des Pfades durchfuehren
*/
          for(i=1;i<=nb;i++) for(j=1;j<=NDIM;j++) x[j][i]+=dx[j][i];
	 } /* Simulationsende */
         erg=erg/nc;
/*
   Querschnitt rausschreiben
*/
	 printf("%f %f %f \n",nb*db,erg,erg*erg);
         fprintf(result,"%f %f %f \n",nb*db,erg,erg*erg);
/*
   r Veteilung rausschreiben
*/
         for(i=0;i<NDIST;i++) 
     fprintf(rdist,"%f %f \n",rmax*(float)i/NDIST,(float)phir[i]/nc);
   /*  fprintf(rdist,"%f %d \n",rmax*(float)i/NDIST,phir[i]);i */
}

float vpot(r)  /* Streupotential */
float r;
{
   return(v0*exp(-r*r/(b*b*2.))) ;
}

float vref(r)  /* Referenzpotential */
float r;
{
   if (r<br)
         return(v0r);
   else
         return(0.);
}


makran() 
{
	extern double ran();

	int i,j;
        float x1,x2;
	float sn,cs,ll;

        for(i=1;i<=nb;i+=2) 
	   for(j=1;j<=NDIM;j++) {
	      do
                x1=ran();
              while (x1<1.0e-7);
	     x2=ran();
             ll=sqrt(-2*log(x1));
             sn=sin(2*pi*x2);
	     cs=cos(2*pi*x2);
             w[j][i]   = ds*ll*cs;
             w[j][i+1] = ds*ll*sn;
          }
}

double ran()
{
    double x;
/*
   THIS ROUTINE CALCULATES A RANDOM NUMBER.
*/
  xk=xk+.5;
  xr=a*xr+floor(c*xk);
  xr=xr-floor(xr/pm)*pm;
  return(xr/pm);
}
