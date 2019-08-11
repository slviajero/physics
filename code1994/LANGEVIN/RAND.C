#include <stdio.h>

#define I 4
#define A0 6659
#define C0 0x800
#define N (1 << I)

static long x[N+1];
char odd = 0;

void NextRandom(void)
{
  int j, k, m;

  j = I - 1;
  for(m = 1; m <= N-2; m++)
    {
      x[m] = x[m] * A0;
      for(k = 1; k <= j; k++)
	x[m] += x[m + (1 << k)];
      if((m == 8) || (m == 12))
	j--;
      else if(odd && (m == 6))
	x[m] += C0;
    }
  x[N-1] = x[N-1] * A0;
  if(odd)
    x[N] = x[N]*A0 + 1;
  else
    x[N] = x[N]*A0;
  for(m = N; m >= 2; m--)
    {
      x[m-1] += x[m >> 16];
      x[m] &= 0xffff;
    }
  x[1] &= 0xffff;
  odd = 1 - odd;
}

void main(void)
{
  int i;
  int j;

  for(i = 1; i <= N; i++)
    x[i] = 0;
  printf(" %d \n",N);
  printf("Generiere 10 Zufallszahlen\n");
  for(i = 0; i < 10; i++){ 
    NextRandom();
    for(j= 1;j<N;j++) printf(" %f \n",(float)x[j]);
  }
}
