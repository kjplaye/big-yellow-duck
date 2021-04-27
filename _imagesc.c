#include <math.h>
#include <SDL.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_LINE 100000
#define MAX_ROWS 100000
#define POWER_STEP 1.2
#define NULL_COLOR 0x305050
#define ZOOM_RATIO 1.2
#define PAN_RATIO 0.3
#define COLOR_MODES 5

#define GET_MAX(a,b) (((a)>(b)) ? (a) : (b))

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
  {1,1,1}
};

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

int imagesc(double * data, long long data_len, long long width, long long frames)
{
  long long i;
  double max, min, amp, power, temp;
  SDL_Event event;
  double win_x0 = 0;
  double win_y0 = 0;
  double win_x1 = 1;
  double win_y1 = 1;
  double xx,yy;
  double alpha,v1,v2;
  long long zz, zz2, zz3;
  long long lines;
  int x,y;
  int r,g,b;
  int flag;
  int redraw_flag = 1;
  int color_mode = 0;
  long long frame_number = 0;
  long long frame_number2;
  long long frame_number3;

  long long last_x,last_y;
  int last_mouse_x = 0, last_mouse_y = 0;
 
  screen_init();

  max = data[0];
  min = data[0];

  for(i=0;i<data_len * frames;i++)
    {
      if (data[i] > max) max = data[i];
      if (data[i] < min) min = data[i];
    }

  power = 1;

  flag = 1;
  redraw_flag = 1;
  while(flag)
    {
      //Redraw
      if (redraw_flag)
	{
	  lines = ceil((double)data_len / width);
	  for(x=0;x<SCREEN_WIDTH;x++)
	    for(y=0;y<SCREEN_HEIGHT;y++)
	      {
		xx = (x * (win_x1 - win_x0)) / SCREEN_WIDTH + win_x0;
		yy = (y * (win_y1 - win_y0)) / SCREEN_HEIGHT + win_y0;
		zz = (long long)(yy * lines) * width + (long long)(xx * width);
		if (xx < 0 | yy < 0 | xx >= 1 | zz >= data_len)
		  point(x,y) = NULL_COLOR;
		else
		  {
		    zz2 = zz;
		    zz += (long long) frame_number * data_len;
		    amp = (data[zz] - min) / (max - min);	   
		    amp = pow(amp,power);
		    switch(color_mode)
		      {
		      case 0:
			r = 255 * amp;
			g = 255 * amp;
			b = 255 * amp;
			break;
		      case 1:
			i = floor(amp * 6);
			v2 = amp*6 - i;
			r = 255 * (pal[i][0]*(1 - v2) + pal[i+1][0]*v2);
			g = 255 * (pal[i][1]*(1 - v2) + pal[i+1][1]*v2);
			b = 255 * (pal[i][2]*(1 - v2) + pal[i+1][2]*v2);
			break;
		      case 2:

			frame_number2 = frame_number + 1;
			if (frame_number2 >= frames) frame_number2 = 0;

			frame_number3 = frame_number2 + 1;
			if (frame_number3 >= frames) frame_number3 = 0;
			
			zz3 = zz2 + (long long) frame_number * data_len;
			amp = (data[zz3] - min) / (max - min);	   
			amp = pow(amp,power);
			r = 255 * amp;

			zz3 = zz2 + (long long) frame_number2 * data_len;
			amp = (data[zz3] - min) / (max - min);	   
			amp = pow(amp,power);
			g = 255 * amp;

			zz3 = zz2 + (long long) frame_number3 * data_len;
			amp = (data[zz3] - min) / (max - min);	   
			amp = pow(amp,power);
			b = 255 * amp;
			break;
		      case 3:
			r = 255 * (1 - amp);
			g = 255 * (1 - amp);
			b = 255 * (1 - amp);
			break;
		      case 4:
			r=g=b=0;
			temp = GET_MAX(max,-min);
			temp *= 1.1;
			if (data[zz] > 0)
			  {
			    amp = (data[zz] - 0) / temp;
			    amp = pow(amp,power);
			    r = 255 * amp;
			  }
			else
			  {
			    amp = (0-data[zz]) / temp;
			    amp = pow(amp,power);
			    b = 255 * amp;
			  }
			break;
		      }
		    point(x,y) = (r << 16) ^ (g << 8) ^ b;
		  }
	      }
	  refresh();
	  redraw_flag = 0;
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
	      redraw_flag = 1;
	      break;
	    case SDL_KEYDOWN:
	      switch(event.key.keysym.sym)
		{
		case SDLK_q:
		  flag = 0;
		  break;
		case SDLK_s:
		  win_x0 = win_y0 = 0;
		  win_x1 = win_y1 = 1;
		  redraw_flag = 1;
		  break;
		case SDLK_c:
		  if (++color_mode >= COLOR_MODES) color_mode = 0;
		  printf("Color mode = %d\n",color_mode);
		  redraw_flag = 1;
		  break;
		case SDLK_PAGEUP:
		  v1 = (win_x0 + win_x1) / 2.0;
		  v2 = (win_x1 - win_x0) * ZOOM_RATIO;
		  win_x0 = v1 - v2*0.5;
		  win_x1 = v1 + v2*0.5;
		  redraw_flag = 1;
		  break;
		case SDLK_PAGEDOWN:
		  v1 = (win_x0 + win_x1) / 2.0;
		  v2 = (win_x1 - win_x0) / ZOOM_RATIO;
		  win_x0 = v1 - v2*0.5;
		  win_x1 = v1 + v2*0.5;
		  redraw_flag = 1;
		  break;
		case SDLK_HOME:
		  v1 = (win_y0 + win_y1) / 2.0;
		  v2 = (win_y1 - win_y0) * ZOOM_RATIO;
		  win_y0 = v1 - v2*0.5;
		  win_y1 = v1 + v2*0.5;
		  redraw_flag = 1;
		  break;
		case SDLK_END:
		  v1 = (win_y0 + win_y1) / 2.0;
		  v2 = (win_y1 - win_y0) / ZOOM_RATIO;
		  win_y0 = v1 - v2*0.5;
		  win_y1 = v1 + v2*0.5;
		  redraw_flag = 1;
		  break;
                case SDLK_KP_PLUS:
                case SDLK_PLUS:
                case SDLK_EQUALS:
                    
		  power /= POWER_STEP;
		  redraw_flag = 1;
		  break;
                case SDLK_KP_MINUS:
                case SDLK_MINUS:
		  power *= POWER_STEP;
		  redraw_flag = 1;
		  break;
		case SDLK_LEFTBRACKET:
		  width++;
		  if (width > data_len) width = data_len;
		  redraw_flag = 1;
		  printf("width = %lld\n",width);
		  break;
		case SDLK_RIGHTBRACKET:
		  width--;
		  if (width < 1) width = 1;
		  printf("width = %lld\n",width);
		  redraw_flag = 1;
		  break;
		case SDLK_LEFT:
		  v2 = (win_x1 - win_x0) * PAN_RATIO;
		  win_x0 -= v2;
		  win_x1 -= v2;
		  redraw_flag = 1;
		  break;
		case SDLK_RIGHT:
		  v2 = (win_x1 - win_x0) * PAN_RATIO;
		  win_x0 += v2;
		  win_x1 += v2;
		  redraw_flag = 1;
		  break;
		case SDLK_UP:
		  v2 = (win_y1 - win_y0) * PAN_RATIO;
		  win_y0 -= v2;
		  win_y1 -= v2;
		  redraw_flag = 1;
		  break;
		case SDLK_DOWN:
		  v2 = (win_y1 - win_y0) * PAN_RATIO;
		  win_y0 += v2;
		  win_y1 += v2;
		  redraw_flag = 1;
		  break;
		case SDLK_RETURN:
		  xx = (last_mouse_x * (win_x1 - win_x0)) / SCREEN_WIDTH + win_x0;
		  yy = (last_mouse_y * (win_y1 - win_y0)) / SCREEN_HEIGHT + win_y0;
		  zz = (long long)(yy * lines) * width + (long long)(xx * width);
		  if (xx < 0 | yy < 0 | xx >= 1 | zz >= data_len)
		    {
		      last_x = -1;
		      last_y = -1;
		    }
		  else
		    {
		      zz += (long long) frame_number * data_len;
		      last_x = (long long)(xx * width);
		      last_y = (long long)(yy * lines);
		    }
		  if (last_x >= 0)
		    {		     
		      printf("data[%lld][%lld] = %f (index %lld) (frame %lld)\n",last_y,last_x,data[zz],zz,frame_number);
		    }
		  else
		    printf("Nothing here\n");
		  break;
		case SDLK_COMMA:
		  if (--frame_number<0) frame_number = 0;		  
		  printf("Frame = %lld\n",frame_number);
		  redraw_flag = 1;
		  break;
		case SDLK_PERIOD:
		  if (++frame_number>=frames) frame_number = frames - 1;
		  printf("Frame = %lld\n",frame_number);
		  redraw_flag = 1;
		  break;
		}
	      break;
	    case SDL_MOUSEMOTION:
	      last_mouse_x = event.motion.x;
	      last_mouse_y = event.motion.y;
	      break;
	    case SDL_MOUSEBUTTONDOWN:
	      switch(event.button.button)
		{
		case SDL_BUTTON_LEFT:
		  alpha = 1/ZOOM_RATIO;
		  break;
		case SDL_BUTTON_RIGHT:
		  alpha = ZOOM_RATIO;
		  break;
		case SDL_BUTTON_MIDDLE:
		  printf("PRESS RETURN...laptop has ackward middle or something\n");
		  break;
		case SDL_BUTTON_WHEELUP:
		  alpha = 1/ZOOM_RATIO;
           	  break;
		case SDL_BUTTON_WHEELDOWN:
		  alpha = ZOOM_RATIO;
           	  break;
		}
	      v1 = (double)event.button.x*(win_x1-win_x0)+SCREEN_WIDTH*win_x0;
	      v2 = (win_x1-win_x0)*alpha;
	      win_x1 = -v1 + ((double)event.button.x-SCREEN_WIDTH)*v2;
	      win_x0 = -v1 + event.button.x*v2;
	      win_x0/=-SCREEN_WIDTH;
	      win_x1/=-SCREEN_WIDTH;
	      
	      
	      v1 = (double)event.button.y*(win_y1-win_y0)+SCREEN_HEIGHT*win_y0;
	      v2 = (win_y1-win_y0)*alpha;
	      win_y1 = -v1 + ((double)event.button.y-SCREEN_HEIGHT)*v2;
	      win_y0 = -v1 + event.button.y*v2;	  
	      win_y0/=-SCREEN_HEIGHT;
	      win_y1/=-SCREEN_HEIGHT;
	      	      
	      redraw_flag = 1;
	      break;
	    case SDL_QUIT:
	      flag = 0;
	      break;	
	    }	  
	}
    }

  SDL_Quit();
}  
