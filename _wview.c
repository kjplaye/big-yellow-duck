#include <math.h>
#include <SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define MAX_FFT_BITS 10
#define MIN_FFT_BITS 4
#define BUFFER_SIZE (1<<MAX_FFT_BITS)
#define TEMP_FILE "/home/kevin/DATA/temp.S16_LE"
#define TWO_PI 6.283185307179586476925286766
#define SCALE_RATIO 1.1
#define ORDER 10
#define COLOR_THRESH 192
#define MAX_ORDER 20

#define MODES 2
#define MODE_TIME 0
#define MODE_SPECT 1
int mode = MODE_TIME;
char mode_string[2][100] = {"Time","Spectrogram"};

#define VOWEL_ITER 8000
#define VOWEL_ANIMATE 100
#define VOWEL_BW_THRESH 300
#define COLOR_SCHEMES 4
#define CS_BLACK_ON_WHITE 0
#define CS_WHITE_ON_BLACK 1
#define CS_HEAT 2
#define CS_PHASE 3
int color_scheme = CS_BLACK_ON_WHITE;

#define HZ2MEL(f) (1127.01048 * log(1 + (f)/700.0))
#define MEL2HZ(m) (700.0 * (exp((m) / 1127.01048) - 1))
#define HZ2INDEX(f) ((f) * (1 << FFT_BITS) / 8000.0)
#define INDEX2HZ(i) ((i) * 8000.0 / (1 << FFT_BITS))
#define FREQ2WFREQ(f) ((warp_flag) ? MEL2HZ(2145.0*(f)/4000.0) : (f))
#define FREQ2WFREQ(f) ((warp_flag) ? MEL2HZ(2145.0*(f)/4000.0) : (f))
#define WFREQ2FREQ(f) ((warp_flag) ? (4000.0*HZ2MEL(f)/2145.0) : (f))
#define FREQ2INDEX(f) ((int)HZ2INDEX(FREQ2WFREQ((f))))


//For pctolsp
#define MAXORD	24
#define N	128
#define NB	15
#define EPS	1.e-6
#define FALSE	0
#define TRUE	1

//For pctoroots
#define MAX_TRY 100
#define MAX_ITER 100
#define ROOT_EPS 0.001

#define point(x,y) pnt[(int)(x)+(int)(y)*SCREEN_WIDTH]
SDL_Surface *screen;
unsigned * pnt;

int SCREEN_WIDTH = 800;
int SCREEN_HEIGHT = 600;

int pal_colors = 6;
double pal[7][3] = {
  {0,0,0},
  {0,0,1},
  {0,1,1},
  {0,1,0},
  {1,1,0},
  {1,0,0},
  {0.7,0.7,1}
};


int pctoroots(double * lpc,complex * root)
{
  //Newton's method: 
  //  x ---> x - f(x)/f'(x)
  int i,j,k;

  int roots_found = 0;
  complex x,last_x,pow_x,f,f_prime;
  double dist;

  for (i=0;i<MAX_TRY && roots_found < 10;i++)
    {
      x = ((4.0*(rand() % RAND_MAX)/RAND_MAX) - 2) + ((4.0*(rand() % RAND_MAX)/RAND_MAX) - 2)*I;
      for(j=0;j<MAX_ITER;j++)
	{
	  last_x = x;

	  pow_x = 1;
	  f = 1;
	  f_prime = 0;
	  for(k=0;k<ORDER;k++)
	    {	      
	      f_prime += pow_x * lpc[k] * (k+1);
	      pow_x *= x;
	      f += pow_x * lpc[k];
	    }
	  x = x - f / f_prime;

	  dist = cabs(x - last_x);

	  if (dist > 1/ROOT_EPS) break;
	  if (dist < ROOT_EPS)
	    {
	      if (cimag(x) < 0) x = creal(x) - cimag(x) * I;
	      for (k=0;k<roots_found;k++)
		{
		  if (cabs(x - root[k]) < ROOT_EPS) break;
		}
	      if (k==roots_found) root[roots_found++] = x;
	      break;
	    }
	}
    }
  return roots_found;
}

void pctolsp2(double * lpc,int m,double * freq,int * lspflag)
{
  static double lastfreq[MAXORD];
  double p[MAXORD], q[MAXORD], pi, ang, fm, tempfreq;
  double fr, pxr, tpxr, tfr, pxm, pxl, fl, qxl, tqxr;
  double qxm, qxr, tqxl;
  int mp, mh, nf, mb, jc, i, j;
  double a[MAXORD+1];

  for(i=0;i<ORDER;i++) a[i+1] = lpc[i];
  a[0] = 1;
  
  pi = atan(1.) * 4.0;
  mp = m + 1;
  mh = m / 2;
  
  /* *generate p and q polynomials				 	*/

  for (i = 0; i < mh; i++)
    {
      p[i] = a[i+1] + a[m-i];
      q[i] = a[i+1] - a[m-i];
    }
  
  /* *compute p at f=0.							*/
  
  fl = 0.;
  for (pxl = 1.0, j = 0; j < mh; j++)  pxl += p[j];
  
  /* *search for zeros of p						*/
  
  nf = 0;
  for (i = 1; i <= N; pxl = tpxr, fl = tfr, i++)
    {
      mb = 0;
      fr = i * (0.5 / N);
      pxr = cos(mp * pi * fr);
      for (j = 0; j < mh; j++)
	{
	  jc = mp - (j+1)*2;
	  ang = jc * pi * fr;
	  pxr += cos(ang) * p[j];
	}
      tpxr = pxr;
      tfr = fr;
      if (pxl * pxr > 0.0) continue;

      do
	{
	  mb++;
	  fm = fl + (fr-fl) / (pxl-pxr) * pxl;
	  pxm = cos(mp * pi * fm);
	  
	  for (j = 0; j < mh; j++)
	    {
	      jc = mp - (j+1) * 2;
	      ang = jc * pi * fm;
	      pxm += cos(ang) * p[j];
	    }
	  (pxm*pxl > 0.0) ? (pxl = pxm, fl = fm) : (pxr = pxm, fr = fm);
	  
	} 
      while ((fabs(pxm) > EPS) && (mb < 4));
      
      if ((pxl-pxr) * pxl == 0) 
	{
	  for (j = 0; j < m; j++) freq[j] = (j+1) * 0.04545;
	  printf("pctolsp2: default lsps used, avoiding /0\n");
	  return;
	}
      freq[nf] = fl + (fr-fl) / (pxl-pxr) * pxl;
      nf += 2;
      if (nf > m-2) break;
    }
  
  
  /* *search for the zeros of q(z)					*/
  
  freq[m] = 0.5;
  fl = freq[0];
  qxl = sin(mp * pi * fl);
  for (j = 0; j < mh; j++)
    {
      jc = mp - (j+1) * 2;
      ang = jc * pi * fl;
      qxl += sin(ang) * q[j];
    }
  
  for (i = 2; i < mp; qxl = tqxr, fl = tfr, i += 2)
    {
      mb = 0;
      fr = freq[i];
      qxr = sin(mp * pi * fr);
      for (j = 0; j < mh; j++)
	{
	  jc = mp - (j+1) * 2;
	  ang = jc * pi * fr;
	  qxr += sin(ang) * q[j];
	}
      tqxl = qxl;
      tfr = fr;
      tqxr = qxr;
      
      do
	{
	  mb++;
	  fm = (fl+fr) * 0.5;
	  qxm = sin(mp * pi * fm);
	  
	  
	  for (j = 0; j < mh; j++)
	    {
	      jc = mp - (j+1) * 2;
	      ang = jc * pi * fm;
	      qxm += sin(ang) * q[j];
	    }
	  (qxm*qxl > 0.0) ? (qxl = qxm, fl = fm) : (qxr = qxm, fr = fm);
	  
	} 
      while ((fabs(qxm) > EPS*tqxl) && (mb < NB));
      
      if ((qxl-qxr) * qxl == 0)
	{
	  for (j = 0; j < m; j++) freq[j] = lastfreq[j];
	  printf("pctolsp2: last lsps used, avoiding /0\n");
	  return;
	}
      freq[i-1] = fl + (fr-fl) / (qxl-qxr) * qxl;
    }
  
  /* *** ill-conditioned cases						*/
  
  *lspflag = FALSE;
  if (freq[0] == 0.0 || freq[0] == 0.5) *lspflag = TRUE;
  for (i = 1; i < m; i++)
    {
      if (freq[i] == 0.0 || freq[i] == 0.5) *lspflag = TRUE;
      
      /* *reorder lsps if non-monotonic					*/
      
      if (freq[i]  <  freq[i-1]) 
	{
	  *lspflag = TRUE;
	  printf("pctolsp2: non-monotonic lsps\n");
	  tempfreq = freq[i];
	  freq[i] = freq[i-1];
	  freq[i-1] = tempfreq;
	}
    }
  
  /* *if non-monotonic after 1st pass, reset to last values		*/
  
  for (i = 1; i < m; i++)
    {
      if (freq[i]  <  freq[i-1])
	{
	  printf("pctolsp2: Reset to previous lsp values\n");
	  for (j = 0; j < m; j++) freq[j] = lastfreq[j];
	  break;
	}
    }
  
  for (i = 0; i < m; i++) lastfreq[i] = freq[i];
}


double window_size = 128;
double fft_scale = 100.0;
int FFT_BITS = 10;

complex A[(1<<MAX_FFT_BITS)];
complex B[(1<<MAX_FFT_BITS)];
complex ROU[(1<<MAX_FFT_BITS)];

complex * fft_temp = A;
complex * fft_vect = B;
// Input in fft_vect, Output in fft_vect
void fft_init(void)
{
  int i;

  for(i=0;i<(1<<FFT_BITS);i++) ROU[i] = cexp(TWO_PI*I*i/(1<<FFT_BITS));
}

void fft(void)
{
  int a,b,c;
  int i,k,pivot;
  int two_to_bits = (1<<FFT_BITS);
  int two_to_bits_mo = (1<<FFT_BITS) - 1;
  int pivot_array[MAX_FFT_BITS];
  complex * temp;

  for(k=0;k<MAX_FFT_BITS;k++) pivot_array[k] = 1<<(FFT_BITS-k-1);

  for(k=0;k<FFT_BITS;k++)
    {
      pivot = pivot_array[k];
      for(i=0;i<two_to_bits;i++)
	{
	  a = i & (pivot - 1);
	  b = i - a;
	  c = (a + (b<<1)) & two_to_bits_mo;

	  fft_temp[i] = fft_vect[c] + fft_vect[c+pivot]*ROU[b];
	}
      temp = fft_temp;
      fft_temp = fft_vect;
      fft_vect = temp;
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
      memcpy(b,temp,sizeof(double) * ORDER);

      eps = 0;
      for(i=1;i<k;i++) eps += lpc[i - 1] * ac[k-i];
      lpc[k-1] = 0;
      for(i=0;i<k;i++) lpc[i] += (ac[k] -  eps) * b[i];
    }

  for(i=0;i<ORDER;i++) lpc[i]*=-1;
}


void window_init(double * window)
{
  int i;

  for(i=0;i<BUFFER_SIZE;i++) 
    {
      window[i] = (i < window_size) ?  0.53836 - 0.46164 * cos((TWO_PI * i)/(window_size-1)) : 0;
    }
}

void refresh(void)
{
  if (SDL_MUSTLOCK(screen)) SDL_UnlockSurface(screen);
  SDL_UpdateRect(screen,0,0,0,0);
  if (SDL_MUSTLOCK(screen)) SDL_LockSurface(screen);
}

void screen_init(void)
{
  SDL_Init(SDL_INIT_VIDEO);
  atexit(SDL_Quit);
  
  screen = SDL_SetVideoMode(SCREEN_WIDTH,SCREEN_HEIGHT, 32, SDL_SWSURFACE | SDL_RESIZABLE);
  pnt = screen->pixels ;
  
  SDL_EnableKeyRepeat(SDL_DEFAULT_REPEAT_DELAY,SDL_DEFAULT_REPEAT_INTERVAL);
  
  if ( SDL_MUSTLOCK(screen) ) SDL_LockSurface(screen);
}

int cmp_complex(complex * a, complex * b)
{
  double a_arg = carg(*a);
  double b_arg = carg(*b);

  if (a_arg < 0) a_arg += TWO_PI;
  if (b_arg < 0) b_arg += TWO_PI;

  if (a_arg > b_arg) return 1;
  if (a_arg < b_arg) return -1;
  return 0;
}


struct filter_t
{
  double coef[MAX_ORDER];
  double memory[MAX_ORDER];
  int order;
};

void filter_init(struct filter_t * filt, int order, double * coef)
{ 
  int i;

  for(i=0;i<=order;i++) 
    {
      filt->coef[i] = coef[i];
      filt->memory[i] = 0;
    }
  filt->order = order;
}

void apply_zero_filter(struct filter_t * filt, double * buffer, int size)
{
  int i,j;
  double temp;

  for (i=0;i<size;i++)
    {
      filt->memory[0] = buffer[i];
      temp = 0;
      for (j=filt->order;j>0;j--)
	{
	  temp += filt->memory[j] * filt->coef[j];
	  filt->memory[j] = filt->memory[j-1];
	}
      buffer[i] = temp + filt->memory[0] * filt->coef[0];
    }
}

void apply_pole_filter(struct filter_t * filt, double * buffer, int size)
{
  int i,j;

  if (filt->coef[0] != 1)
    {
      printf("apply_pole_filter:  bad coefficients");
      exit(1);
    }

  for(i=0;i<size;i++)
    {
      filt->memory[0] = buffer[i];
      for(j=filt->order;j>0;j--)
	{
	  filt->memory[0] -= filt->coef[j] * filt->memory[j];
	  filt->memory[j]  = filt->memory[j-1];
	}
      buffer[i] = filt->memory[0];
    }
}


void wview(double * buffer, long size, int frames)
{
  SDL_Event event;
  long i,j,k,flag,x,y,old_y,r,g,b,c;
  long start,end,new_start,new_end;

  int toggle = 1;
  int line_flag = 0;
  int root_flag = 0;
  int grid_flag = 0;
  int hpf_flag = 1;
  int warp_flag = 0;
  
  int selection_event = 1;
  int redraw_event = 1;

  char cmd_line[1000];
  double window[BUFFER_SIZE];
  double small_buffer[BUFFER_SIZE];
  double amp,mean,var,amp_max,phase,phase_offset;
  double min[frames];
  double max[frames];

  double ac[ORDER + 1];
  double lpc[ORDER];

  double amp_1[SCREEN_HEIGHT];
  double amp_2[SCREEN_HEIGHT];

  int lspflag;
  double lsp[ORDER];
  complex root[ORDER];
  int roots_found;
  double formant[2];
  int frame = 0;

  // 2nd order Butterworth 100 Hz HPF:
  struct filter_t hpf_zero;
  struct filter_t hpf_pole;
  double ahpf[3] = {0.946, -1.892, 0.946};
  double bhpf[3] = {1.0, -1.889033, 0.8948743};

  if (size <= 0) return;
  if (frames <= 0) return;

  screen_init();
  refresh();
  fft_init();
  window_init(window);

  start = 0;
  end = size;

  new_start = 0;
  new_end = size-1;

  for(j=0;j<frames;j++)
    {
      min[j] = max[j] = buffer[size*j + 0];
      for(i=0;i<size;i++)
	{
	  if (buffer[size*j + i] > max[j]) max[j] = buffer[size*j + i];
	  if (buffer[size*j + i] < min[j]) min[j] = buffer[size*j + i];      
	}
    }
  
  flag = 1;
  while(flag)
    {        
      //Refresh bounds
      if (start > end)
	{
	  start ^= end;
	  end ^= start;
	  start ^= end;
	}

      if (new_start > new_end)
	{
	  new_start ^= new_end;
	  new_end ^= new_start;
	  new_start ^= new_end;
	}

      if (selection_event)
	{
	  printf("Showing samples %ld..%ld selected %ld..%ld\n",start,end,new_start,new_end);
	  selection_event = 0;
	}

      if (redraw_event)
	{
	  //Clear screen
	  for(x=0;x<SCREEN_WIDTH;x++)
	    for(y=0;y<SCREEN_HEIGHT;y++) point(x,y)=0;
	  
	  //Begin Drawing
	  switch(mode)
	    {
	    case MODE_TIME:	 	      
	      for (i=start;i<end;i++) 
		{
		  x = ((i-start)*SCREEN_WIDTH)/(end-start);
		  y = (1 - ((buffer[size * frame + i] - min[frame]) / (max[frame]-min[frame]))) * SCREEN_HEIGHT;
		  if (y<0) y=0;
		  if (y>=SCREEN_HEIGHT) y = SCREEN_HEIGHT-1;
		  point(x,y) = 0xffffff;
		  if (line_flag && i>start)
		    {
		      if (old_y < y)
			{
			  for(k=old_y;k<y;k++) point(x,k) = 0x808080;
			}
		      if (old_y > y)
			{
			  for(k=y;k<old_y;k++) point(x,k) = 0x808080;
			}		  
		    }
		  old_y = y;
		}
	      break;
	    case MODE_SPECT:
	      for(x=0;x<SCREEN_WIDTH;x++)
		{
		  i = x * (end - start) / SCREEN_WIDTH + start;
		  for(j=0;j<BUFFER_SIZE;j++) 
		    {
		      if (i+j >= size)
			small_buffer[j] = 0;
		      else
			small_buffer[j] = buffer[size * frame + i+j];
		    }
		  phase_offset = i+(window_size - 1)*0.5;		 
		  if (hpf_flag) 
		    {
		      filter_init(&hpf_zero,2,ahpf);
		      filter_init(&hpf_pole,2,bhpf);
		      apply_zero_filter(&hpf_zero,small_buffer,BUFFER_SIZE);
		      apply_pole_filter(&hpf_pole,small_buffer,BUFFER_SIZE);
		    }
		  for(j=0;j<BUFFER_SIZE;j++) fft_vect[j] = small_buffer[j] * window[j];
		  switch (mode)
		    {
		    case MODE_SPECT:
		      fft();
		      break;
		    }
		  
		  for(y=0;y<SCREEN_HEIGHT;y++)
		    {
		      j = FREQ2INDEX(y*4000.0/SCREEN_HEIGHT);			
		      amp = fft_scale * cabs(fft_vect[j]);
		      i = x * (end - start) / SCREEN_WIDTH + start; 
		      phase = carg(fft_vect[j] * cexp(-TWO_PI*I*phase_offset/(1 << (FFT_BITS-1))));
		      
		      if (amp >= COLOR_THRESH) amp = 255 - (255*COLOR_THRESH - COLOR_THRESH*COLOR_THRESH)/amp; 
		      switch(color_scheme)
			{
			case CS_BLACK_ON_WHITE:
			  i = amp;
			  if (i >= 256) i = 255;
			  point(x,SCREEN_HEIGHT - 1 - y) = (255 - i) * 0x010101;
			  break;
			case CS_WHITE_ON_BLACK:
			  i = amp;
			  if (i >= 256) i = 255;
			  point(x,SCREEN_HEIGHT - 1 - y) = (i) * 0x010101;
			  break;
			case CS_HEAT:
			  amp = amp*pal_colors/256.0;
			  i = amp;
			  if (i>=pal_colors) {i=pal_colors-1;amp=pal_colors;}	 

			  r = 256*(pal[i+1][0]*(amp-i) + pal[i][0]*(i+1-amp));
			  g = 256*(pal[i+1][1]*(amp-i) + pal[i][1]*(i+1-amp));
			  b = 256*(pal[i+1][2]*(amp-i) + pal[i][2]*(i+1-amp));
			  if (r>255) r = 255;
			  if (g>255) g = 255;
			  if (b>255) b = 255;

			  if (r<0) r=0;
			  if (g<0) g=0;
			  if (b<0) b=0;
			  point(x,SCREEN_HEIGHT - 1 - y) = r*0x10000 + g*0x100 + b;
			  break;
			case CS_PHASE:
			  r = g = b = 1.15470054;

			  r+=cos(phase) * -0.70710678119;
			  g+=cos(phase) *  0.70710678119;

			  r+=sin(phase) *  0.40824829046;
			  g+=sin(phase) *  0.40824829046;
			  b+=sin(phase) * -0.81649658093;		  
			  
			  r*=amp;
			  g*=amp;
			  b*=amp;
			  
			  if (r<0) r=0;
			  if (g<0) g=0;
			  if (b<0) b=0;

			  if (r>255) r = 255;
			  if (g>255) g = 255;
			  if (b>255) b = 255;

			  point(x,SCREEN_HEIGHT - 1 - y) = (r<<16) + (g<<8) + b;
			  break;
			}
		    }
		}
	      break;
	    }
	  //Draw bounds
	  x = ((new_start-start)*SCREEN_WIDTH)/(end-start);      
	  if (x>=0 && x<SCREEN_WIDTH)
	    {
	      for(y=0;y<SCREEN_HEIGHT;y++) point(x,y)^=0x00ff00;
	    }
	  
	  x = ((new_end-start)*SCREEN_WIDTH)/(end-start);      
	  if (x>=0 && x<SCREEN_WIDTH)
	    {
	      for(y=0;y<SCREEN_HEIGHT;y++) point(x,y)^=0xff0000;
	    }

	  if (grid_flag && mode!=MODE_TIME)
	    {
	      for(x=0;x<SCREEN_WIDTH;x++) 
		{
		  point(x,(SCREEN_HEIGHT * 1) / 4)^=0x404040;
		  point(x,(SCREEN_HEIGHT * 2) / 4)^=0x404040;
		  point(x,(SCREEN_HEIGHT * 3) / 4)^=0x404040;
		}
	    }
	  
	  refresh();
	  redraw_event = 0;
	}

      //Event loop
      if (SDL_WaitEvent(&event))
	{
	  switch(event.type)
	    {
	    case SDL_VIDEORESIZE:
	      SCREEN_WIDTH = event.resize.w;
	      SCREEN_HEIGHT = event.resize.h;
	      screen_init();
	      redraw_event = 1;
	      break;
	    case SDL_KEYDOWN:
	      if (event.key.keysym.sym == SDLK_LEFTBRACKET) 
		{

		  x = (end - start)/3;
		  start -= x;
		  end -= x;
		  new_start -= x;
		  new_end -= x;

		  if (start < 0) 
		    {
		      end+=-start;
		      new_end+=-start;
		      new_start+=-start;
		      start+=-start;
		    }
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_RIGHTBRACKET) 
		{
		  x = (end - start)/3;
		  start += x;
		  end += x;
		  new_start += x;
		  new_end += x;

		  if (end > size) 
		    {
		      start -= end - size;
		      new_end -= end - size;
		      new_start -= end - size;
		      end -= end - size;
		    }
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_UP) 
		{
		  fft_scale *= SCALE_RATIO;
		  printf("new fft scale of %f\n",fft_scale);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_DOWN) 
		{
		  fft_scale /= SCALE_RATIO;
		  printf("new fft scale of %f\n",fft_scale);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_RIGHT) 
		{
		  window_size *= SCALE_RATIO;
		  if (window_size > (1<<FFT_BITS)) window_size = (1<<FFT_BITS);
		  printf("new window size of %f\n",window_size);
		  window_init(window);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_LEFT) 
		{
		  window_size /= SCALE_RATIO;
		  printf("new window size of %f\n",window_size);
		  window_init(window);
		  redraw_event = 1;
		}

	      if (event.key.keysym.sym == SDLK_PAGEUP) 
		{
		  if (++FFT_BITS > MAX_FFT_BITS) FFT_BITS = MAX_FFT_BITS;
		  fft_init();
		  printf("fft is now %d bits\n",FFT_BITS);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_PAGEDOWN) 
		{
		  if (--FFT_BITS < MIN_FFT_BITS) FFT_BITS = MIN_FFT_BITS;
		  fft_init();
		  printf("fft is now %d bits\n",FFT_BITS);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_m) 
		{
		  if (++mode == MODES) mode=0;
		  printf("Entering mode: %s\n",mode_string[mode]);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_c) 
		{
		  if (++color_scheme == COLOR_SCHEMES) color_scheme = 0;
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_v)
		{
		  for(x=0;x<SCREEN_WIDTH;x++)
		    for(y=0;y<SCREEN_HEIGHT;y++) point(x,y) = 0;

		  for(x=0;x<SCREEN_WIDTH;x++)
		    {
		      point(x,1*SCREEN_HEIGHT/(4)) = 0x0000ff;
		      point(x,2*SCREEN_HEIGHT/(4)) = 0x0000ff;
		      point(x,3*SCREEN_HEIGHT/(4)) = 0x0000ff;
		    }
		  for(y=0;y<SCREEN_HEIGHT;y++)
		    {
		      point(1*SCREEN_WIDTH/(4),y) = 0x0000ff;
		      point(2*SCREEN_WIDTH/(4),y) = 0x0000ff;
		      point(3*SCREEN_WIDTH/(4),y) = 0x0000ff;
		    }

		  c = 0;
		  for(i=new_start;i<new_end;i+=(1+(new_end - new_start)/VOWEL_ITER))
		    {
 		      for(j=0;j<BUFFER_SIZE;j++) 
			{
			  if (i+j >= size)
			    small_buffer[j] = 0;
			  else
			    small_buffer[j] = buffer[size * frame + i+j];
			}
		      if (hpf_flag)
			{
			  filter_init(&hpf_zero,2,ahpf);
			  filter_init(&hpf_pole,2,bhpf);
			  apply_zero_filter(&hpf_zero,small_buffer,BUFFER_SIZE);
			  apply_pole_filter(&hpf_pole,small_buffer,BUFFER_SIZE);
			}
		      for(j=0;j<BUFFER_SIZE;j++) fft_vect[j] = small_buffer[j] * window[j];

		      //Find auto correlations
		      for(k=0;k<=ORDER;k++)
			{
			  ac[k] = 0;
			  for(j=0;j<BUFFER_SIZE - k;j++) ac[k] += creal(fft_vect[j])*creal(fft_vect[j+k]);
			}
		      
		      //Do Durbin recursion
		      durbin(ac,lpc);
		      roots_found = pctoroots(lpc,root);
		      qsort(root,10,sizeof(complex),cmp_complex);

		      j = 0;
		      for(k=0;k<10 && j<2;k++)
			{
			  if (carg(root[k]) * 8000.0 / TWO_PI < 200) continue;
			  if (carg(root[k]) * 8000.0 / TWO_PI < 400 && j>0) continue;
			  if (carg(root[k]) * 8000.0 / TWO_PI > 1000 && j==0) continue;
			  if (carg(root[k]) * 8000.0 / TWO_PI > 3000) continue;

			  if ( (2 * (cabs(root[k]) - 1) * 8000 / TWO_PI) > VOWEL_BW_THRESH) continue;
			  formant[j++] = carg(root[k]) * 8000.0 / TWO_PI;			  
			}
		      if (k!=10)
			{
			  x = formant[0] * SCREEN_WIDTH / (1000.0 );
			  y = (4000 - formant[1]) * SCREEN_HEIGHT / (4000.0 );
			  
			  amp = 255 * (2 * (cabs(root[k]) - 1) * 8000 / TWO_PI) / VOWEL_BW_THRESH;
			  
			  point(x,y) = 0x010101 * (int)amp;		     
			}
		      if (!((++c) % VOWEL_ANIMATE)) refresh();
		    }
		  refresh();
		}
	      if (event.key.keysym.sym == SDLK_f) 
		{
		  hpf_flag = !hpf_flag;
		  printf("high pass filter flag is now %d\n",hpf_flag);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_l) 
		{
		  line_flag = !line_flag;
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_w) 
		{
		  warp_flag = !warp_flag;		  
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_g) 
		{
		  grid_flag = !grid_flag;
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_r) 
		{
		  root_flag = !root_flag;
		  redraw_event = 1;
		}
 	      if (event.key.keysym.sym == SDLK_z) 
		{		 
		  start = new_start;
		  end = new_end;
		  selection_event = 1;
		  redraw_event = 1;				  
		}
	      if (event.key.keysym.sym == SDLK_o) 
		{
		  amp = (SCALE_RATIO - 1) * (end - start) / 2.0;
		  start -= amp;
		  end += amp;
		  if (start < 0) start = 0;
		  if (end > size) end = size;		 
		  selection_event = 1;
		  redraw_event = 1;		  
		}
	      if (event.key.keysym.sym == SDLK_i) 
		{
		  amp = (1.0/SCALE_RATIO - 1) * (end - start) / 2.0;
		  start -= amp;
		  end += amp;
		  if (end <= start) end = start + 1;
		  selection_event = 1;
		  redraw_event = 1;				  
		}
	      if (event.key.keysym.sym == SDLK_s) 
		{
		  start = 0;
		  end = size;
		  selection_event = 1;
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_COMMA)
		{
		  if (--frame<0) frame = 0;
		  printf("Frame = %d\n",frame);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_PERIOD)
		{
		  if (++frame>=frames) frame = frames - 1;
		  printf("Frame = %d\n",frame);
		  redraw_event = 1;
		}
	      if (event.key.keysym.sym == SDLK_ESCAPE) flag = 0;
	      if (event.key.keysym.sym == SDLK_TAB) 
		{
		  toggle = !toggle;
		  SDL_WM_ToggleFullScreen(screen);
		}
	      break;
	    case SDL_MOUSEBUTTONDOWN:
	      switch(event.button.button)
		{
		case SDL_BUTTON_LEFT:
		  new_start = event.button.x * (end - start) / SCREEN_WIDTH + start;
		  selection_event = 1;
		  redraw_event = 1;
		  break;
		case SDL_BUTTON_RIGHT:
		  new_end = event.button.x * (end - start) / SCREEN_WIDTH + start;
		  selection_event = 1;
		  redraw_event = 1;				  
		  break;		  
		}
	      break;
	    case SDL_QUIT:
	      flag = 0;
	      break;	
	    }	  
	}
    }
  
  SDL_Quit();
}  
