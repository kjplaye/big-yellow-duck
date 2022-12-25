#include <SDL2/SDL.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <pthread.h>

#define MAX_FFT_BITS 10
#define MIN_FFT_BITS 4
#define BUFFER_SIZE (1<<MAX_FFT_BITS)
#define TEMP_FILE "/home/kevin/DATA/temp.S16_LE"
#define TWO_PI 6.283185307179586476925286766
#define SCALE_RATIO 1.1
#define ORDER 10
#define COLOR_THRESH 192
#define MAX_ORDER 20
#define AUDIO_BUFFER_SIZE 128
#define AUDIO_FOLLOW_RATIO 0.618
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MAX_THREADS 64

#define MODES 2
#define MODE_TIME 0
#define MODE_SPECT 1
int mode = MODE_TIME;
char mode_string[2][100] = {"Time","Spectrogram"};

#define VOWEL_ITER 8000
#define VOWEL_ANIMATE 100
#define VOWEL_BW_THRESH 300
#define COLOR_SCHEMES 3
#define CS_BLACK_ON_WHITE 0
#define CS_WHITE_ON_BLACK 1
#define CS_HEAT 2
int color_scheme = CS_BLACK_ON_WHITE;

#define HZ2MEL(f) (1127.01048 * log(1 + (f)/700.0))
#define MEL2HZ(m) (700.0 * (exp((m) / 1127.01048) - 1))
#define HZ2INDEX(f) ((f) * (1 << FFT_BITS) / 8000.0)
#define INDEX2HZ(i) ((i) * 8000.0 / (1 << FFT_BITS))
#define FREQ2INDEX(f) ((int)HZ2INDEX(f))

#define point(x,y) pnt[(int)(x)+(int)(y)*SCREEN_WIDTH]
SDL_Window *screen;
SDL_Renderer * renderer;
SDL_Texture * texture;
unsigned * pnt;

#define MAX_SCREEN_WIDTH 10000
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

double window_size = 128;
double fft_scale = 100.0;
int FFT_BITS = 10;

complex A[MAX_SCREEN_WIDTH][(1<<MAX_FFT_BITS)];
complex B[MAX_SCREEN_WIDTH][(1<<MAX_FFT_BITS)];
complex ROU[(1<<MAX_FFT_BITS)];

complex * fft_temp[MAX_SCREEN_WIDTH];
complex * fft_vect[MAX_SCREEN_WIDTH];
// Input in fft_vect, Output in fft_vect
void fft_init(void)
{
  int i,x;

  for(i=0;i<(1<<FFT_BITS);i++) ROU[i] = cexp(TWO_PI*I*i/(1<<FFT_BITS));
  for(x=0;x<MAX_SCREEN_WIDTH;x++) fft_temp[x] = A[x];
  for(x=0;x<MAX_SCREEN_WIDTH;x++) fft_vect[x] = B[x];
}

void fft(int x)
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

	  fft_temp[x][i] = fft_vect[x][c] + fft_vect[x][c+pivot]*ROU[b];
	}
      temp = fft_temp[x];
      fft_temp[x] = fft_vect[x];
      fft_vect[x] = temp;
    }
}

void window_init(double * window)
{
  int i;

  for(i=0;i<BUFFER_SIZE;i++) 
    {
      window[i] = (i < window_size) ?  0.53836 - 0.46164 * cos((TWO_PI * i)/(window_size-1)) : 0;
    }
}

void refresh()
{
  SDL_UpdateTexture(texture, NULL, pnt, SCREEN_WIDTH * sizeof(unsigned));
  SDL_RenderClear(renderer);
  SDL_RenderCopy(renderer, texture, NULL, NULL);
  SDL_RenderPresent(renderer);
}

void screen_init(int init)
{
  if (init)
    {
      SDL_Init(SDL_INIT_VIDEO);
      atexit(SDL_Quit);
      screen = SDL_CreateWindow("imagesc",SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
				SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_RESIZABLE);
      renderer = SDL_CreateRenderer(screen, -1, 0);
    }
  
  SDL_RenderClear(renderer);
  SDL_GL_SwapWindow(screen);    
  texture = SDL_CreateTexture(renderer,SDL_PIXELFORMAT_ARGB8888,
			      SDL_TEXTUREACCESS_STREAMING,
			      SCREEN_WIDTH, SCREEN_HEIGHT);
  SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY,
	      "linear");  // make the scaled rendering look smoother.
  SDL_RenderSetLogicalSize(renderer, SCREEN_WIDTH, SCREEN_HEIGHT);

  if (pnt) free(pnt);
  if ((pnt = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(unsigned))) == 0)
    {fprintf(stderr, "OUT OF MEMORY");exit(1);}
}

double * audio_buffer;
int audio_current;
int audio_start;
int audio_stop;
double audio_abs_max;
int audio_window_start;
int audio_window_end;
int audio_window_bar;
int audio_follow;
int audio_playing;
pthread_mutex_t audio_mutex;
void audio_callback(void * userdata, unsigned char * stream, int len)
{
  int i,y;
  unsigned short * data_stream = (unsigned short *)stream;
  int stream_size = len / 2;
  unsigned short sample;
  
  if (audio_current >= audio_stop)
    {
      if (audio_window_bar >= 0 && audio_window_bar < SCREEN_WIDTH)
      	{
      	  for(y=0;y<SCREEN_HEIGHT;y++) point(audio_window_bar,y)^=0xffffff;
      	}
      pthread_mutex_lock(&audio_mutex);
      refresh();
      pthread_mutex_unlock(&audio_mutex);
      audio_follow = 0;
      audio_playing = 0;
      SDL_PauseAudio(1);
      return;
    }

  if (audio_window_bar >= 0 && audio_window_bar < SCREEN_WIDTH)
    {
      for(y=0;y<SCREEN_HEIGHT;y++) point(audio_window_bar,y)^=0xffffff;
    }
  pthread_mutex_lock(&audio_mutex);
  refresh();
  pthread_mutex_unlock(&audio_mutex);
  audio_window_bar = ((audio_current-audio_window_start)*SCREEN_WIDTH)/(audio_window_end-audio_window_start);
  if (audio_window_bar >= 0 && audio_window_bar < SCREEN_WIDTH)
    {
      for(y=0;y<SCREEN_HEIGHT;y++) point(audio_window_bar,y)^=0xffffff;
    }
  pthread_mutex_lock(&audio_mutex);
  refresh();
  pthread_mutex_unlock(&audio_mutex);

  
  for(i=0;i<stream_size;i++)
    {
      sample = 32767 * audio_buffer[audio_current] / audio_abs_max;
      data_stream[i] = (audio_current >= audio_stop) ? 0 : sample;
      audio_current++;
    }
}

void play_audio(double * data, int start, int stop)  
{
  int y;

  if (audio_playing && audio_follow == 0)
    {
      if (audio_window_bar >= 0 && audio_window_bar < SCREEN_WIDTH)
      	{
      	  for(y=0;y<SCREEN_HEIGHT;y++) point(audio_window_bar,y)^=0xffffff;
      	}
      pthread_mutex_lock(&audio_mutex);
      refresh();
      pthread_mutex_unlock(&audio_mutex);

      audio_playing = 0;
      audio_follow = 0;
      SDL_PauseAudio(1);
      return;
    }

  audio_current = start;
  audio_start = start;
  audio_stop = stop;
  audio_window_bar = -1;
  audio_playing = 1;
  SDL_PauseAudio(0);
}

void audio_init(double * data, double min, double max)
{
  SDL_AudioSpec *desired, *obtained;
  SDL_AudioSpec *hardware_spec;

  audio_buffer = data;
  audio_abs_max = MAX(max,-min);

  desired = malloc(sizeof(SDL_AudioSpec));
  obtained = malloc(sizeof(SDL_AudioSpec));
  desired->freq=4000;
  desired->format=AUDIO_U16LSB;
  desired->channels=0;  
  desired->samples=AUDIO_BUFFER_SIZE;
  desired->callback=audio_callback;
  desired->userdata=NULL;
  if ( SDL_OpenAudio(desired, obtained) < 0 ){
    fprintf(stderr, "Couldn't open audio: %s\n", SDL_GetError());
    exit(-1);
  }
  free(desired);
  hardware_spec=obtained;
}

int cmp_complex(const void * aa, const void  * bb)
{
  const complex * a = aa;
  const complex * b = bb;
  
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

struct strip_t
{
  int x;
  long start;
  long end;
  long size;
  int frame;
  double * buffer;
  double * window;
};

void * draw_fft_strip(void * arg)
{
  double small_buffer[BUFFER_SIZE];
  long i,j,y,r,g,b;
  double amp;

  struct strip_t * param = arg;
  
  int x = param->x;
  long start = param->start;
  long end = param->end;
  long size = param->size;
  int frame = param->frame;
  double * buffer = param->buffer;
  double * window = param->window;

  i = x * (end - start) / SCREEN_WIDTH + start;
  for(j=0;j<BUFFER_SIZE;j++) 
    {
      if (i+j >= size)
	small_buffer[j] = 0;
      else
	small_buffer[j] = buffer[size * frame + i+j];
    }
  for(j=0;j<BUFFER_SIZE;j++) fft_vect[x][j] = small_buffer[j] * window[j];
  fft(x);
  
  for(y=0;y<SCREEN_HEIGHT;y++)
    {
      j = FREQ2INDEX(y*4000.0/SCREEN_HEIGHT);			
      amp = fft_scale * cabs(fft_vect[x][j]);
      i = x * (end - start) / SCREEN_WIDTH + start; 
      
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
	}
    }
  return NULL;
}

void wview(double * buffer, long size, int frames)
{
  SDL_Event event;
  long i,j,k,flag,x,y,old_y,r,g,b,c;
  long start,end,new_start,new_end;
  struct strip_t param[MAX_SCREEN_WIDTH];

  int toggle = 1;
  int line_flag = 0;
  int root_flag = 0;
  int grid_flag = 0;
  
  int selection_event = 1;
  int redraw_event = 1;

  char cmd_line[1000];
  double window[BUFFER_SIZE];
  double amp,mean,var,amp_max;
  double min[frames];
  double max[frames];
  pthread_t thread[MAX_THREADS];

  double amp_1[SCREEN_HEIGHT];
  double amp_2[SCREEN_HEIGHT];

  int lspflag;
  double lsp[ORDER];
  complex root[ORDER];
  int roots_found;
  double formant[2];
  int frame = 0;
  double nstart, nend;

  if (size <= 0) return;
  if (frames <= 0) return;

  screen_init(1);
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

  audio_init(buffer, min[0], max[0]);
  
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

      audio_window_start = start;
      audio_window_end = end;
      if (selection_event)
	{
	  printf("Showing samples %ld..%ld selected %ld..%ld\n",start,end,new_start,new_end);
	  selection_event = 0;
	}

      if (redraw_event)
	{
	  // Draw audio bar
	  if (audio_window_bar >= 0 && audio_window_bar < SCREEN_WIDTH)
	    {
	      for(y=0;y<SCREEN_HEIGHT;y++) point(audio_window_bar,y)^=0xffffff;
	    }	  
	  
	  //Begin Drawing
	  switch(mode)
	    {
	    case MODE_TIME:
	      //Clear screen
	      for(x=0;x<SCREEN_WIDTH;x++)
		for(y=0;y<SCREEN_HEIGHT;y++) point(x,y)=0;
	      
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
		  param[x].x = x;
		  param[x].start = start;
		  param[x].end = end;
		  param[x].size = size;
		  param[x].frame = frame;
		  param[x].buffer = buffer;
		  param[x].window = window;

		  pthread_create(&thread[x % MAX_THREADS], NULL, draw_fft_strip, &param[x]);
		  if ((x % MAX_THREADS) == MAX_THREADS - 1)
		    {
		      for(i=0;i<MAX_THREADS;i++) pthread_join(thread[i], NULL);
		    }
		}
	      if (x % MAX_THREADS)
		for(i=0;i<x % MAX_THREADS;i++) pthread_join(thread[i], NULL);
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

	  pthread_mutex_lock(&audio_mutex);	   
	  refresh();
	  pthread_mutex_unlock(&audio_mutex);

	  redraw_event = 0;
	}

      //Event loop
      if (SDL_PollEvent(&event) || audio_follow)
	{
	  if (audio_follow)	   
	    {
	      if ((audio_current - start) / (end - start) > AUDIO_FOLLOW_RATIO)
		{
		  nstart = audio_current - AUDIO_FOLLOW_RATIO * (end - start);
		  nend = nstart + (end - start);
		  start = nstart;
		  end = nend;
		  if (end >= size)
		    {
		      end = size;
		      start = end - (nend - nstart);
		    }
		  selection_event = 1;
		  redraw_event = 1;
		}
	    }	  
	  switch(event.type)
	    {
	    case SDL_WINDOWEVENT:
	      if (event.window.event == SDL_WINDOWEVENT_RESIZED)
		{
		  SCREEN_WIDTH = event.window.data1;
		  SCREEN_HEIGHT = event.window.data2;
		  screen_init(0);
		  redraw_event = 1;
		}
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
	      if (event.key.keysym.sym == SDLK_p) 
		{
		  audio_follow = 0;
		  play_audio(buffer, new_start, new_end);
		}
	      if (event.key.keysym.sym == SDLK_f) 
		{
		  audio_follow = 1;
		  play_audio(buffer, start, size);
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
	      if (event.key.keysym.sym == SDLK_l) 
		{
		  line_flag = !line_flag;
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
