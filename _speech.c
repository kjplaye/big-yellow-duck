#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ORDER 10
#define PULSE_SHAPE 1.3


double randn(void)
{
  static int next_available = 0;
  static double x;
  static double y;
  double t;
  double s;
  
  if (next_available)
    {
      next_available = 0;
      return y;
    }

  while(1)
    {
      x = rand() * 2.0 / RAND_MAX - 1.0;
      y = rand() * 2.0 / RAND_MAX - 1.0;
      t = x*x + y*y;
      if (t<=1) break;
    }
  
  s = sqrt(-2.0 * log(t) / t);
  x *= s;
  y *= s;
 
  next_available = 1;
  return x;
}

void zero_filt_mul(double * mem, double * coef, double *buff, int coef_len, int buff_len)
{
  int i,j;
  double temp;
  for(i=0;i<buff_len;i++) 
    {
      mem[0] = buff[i];
      temp = 0;
      for(j=coef_len-1;j>0;j--)
	{
	  temp += mem[j] * coef[j];
	  mem[j] = mem[j-1];
	}
      buff[i] = temp + mem[0] * coef[0];             
    }
}


void pole_filt_mul(double * mem, double * coef, double *buff, int coef_len, int buff_len)
{
  int i,j;
  for(i=0;i<buff_len;i++)
    {
      mem[0] = buff[i];
      for(j=coef_len-1;j>0;j--)
	{
	  mem[0] -= coef[j] * mem[j];
	  mem[j] = mem[j-1];
	}
      buff[i] = mem[0];
    }
}

void durbin(double * ac, double * lpc)
{
  int i,k;
  double alpha,beta,eps;
  double b[ORDER];
  double temp[ORDER];
  
  if (ac[0] == 0)
    {
      for(i=0;i<ORDER;i++) lpc[i] = 0;
      return;
    }

  lpc[0] = ac[1] / ac[0];
  b[0] = 1 / ac[0];

  for(k=2;k<=ORDER;k++)
    {
      eps = 0;
      for(i=1;i<k;i++) eps += ac[i] * b[i-1];
      alpha = 1/(1 - eps*eps);
      beta = - eps * alpha;

      for(i=k-1;i>0;i--) b[i] = b[i-1];
      b[0] = 0;
      for(i=0;i<k;i++) temp[i] = alpha *  b[i] + beta * b[k - 1 - i];
      for(i=0;i<ORDER;i++) b[i]=temp[i];

      eps = 0;
      for(i=1;i<k;i++) eps += lpc[i - 1] * ac[k-i];
      lpc[k-1] = 0;
      for(i=0;i<k;i++) lpc[i] += (ac[k] -  eps) * b[i];
    }
}

void lpc2cepst(double * lpc, double * cepst)
{
  double temp_cepst[10];
  double temp_lpc[10];
  int i,k;

  for(i=0;i<10;i++)
    {
      cepst[i] = lpc[i];
      for(k=1;k<=i;k++)
	cepst[i] += (k/(i+1.0)) * cepst[k-1] * lpc[i-k];
    }
}

double pulse(double x, double pulse_energy)
{
  return exp(-PULSE_SHAPE*x)*x / pulse_energy;
}

void synth(double * ans, int size, double * tslp,double gain, double delay, double gain_on_pitch,double pulse_energy)
{
  int i;
  double x;

  for(i=0;i<size;i++)
    {
      x = randn();
      *tslp+=1;
      if (*tslp >= delay) *tslp = *tslp - ((int)(*tslp/delay))*delay;
      x += pulse(*tslp,pulse_energy) * gain_on_pitch;
      ans[i] = x;
    }
}
