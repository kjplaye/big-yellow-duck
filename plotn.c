#include <math.h>
#include <SDL.h>
#include <stdio.h>
#include <stdlib.h>

#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 600
#define SCALE_RATIO 1.1
#define MAX_LINE 500000
#define MAX_ROWS 500000
#define MAX_DIM 100
#define POWER_STEP 1.1
#define INDEX(i,j) ((i)*(n) + (j))
#define TWO_PI 6.283185307179586476925286766

#define BACK_COLOR 0x000000
#define PEN_COLOR 0xffffff
#define MAX_COLOR 7

#define point(x,y) pnt[(int)(x)+(int)(y)*SCREEN_WIDTH]
SDL_Surface *screen;
unsigned * pnt;

unsigned pal[MAX_COLOR] = {
  0x8080ff,
  0x80ff80,
  0x80ffff,
  0xff8080,
  0xff80ff,
  0xffff80,
  0xffffff
};

double X_SCALE = 400;
double Y_SCALE = 400;

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

int cols(FILE * fp)
{
  int i;
  int number;
  int number_of_cols;
  int space_flag;
  char line[MAX_LINE];

  number = fread(line,1,MAX_LINE,fp);
  number_of_cols = 0;
  space_flag = 1;
  for(i=0;i<number;i++)
    {
      if (line[i] == 0x0a) break;
      if (isspace(line[i])) 
	space_flag=1;
      else
	{
	  if (space_flag)
	    {
	      space_flag = 0;
	      number_of_cols++;
	    }
	}
    }
  fseek(fp,0,SEEK_SET);

  return number_of_cols;
}

double * data;
unsigned * col;
double max[MAX_DIM];
double min[MAX_DIM];
int n,m;
double A[10][10] = {
  {1,0,0,0,0,0,0,0,0,0},
  {0,1,0,0,0,0,0,0,0,0},
  {0,0,1,0,0,0,0,0,0,0},
  {0,0,0,1,0,0,0,0,0,0},
  {0,0,0,0,1,0,0,0,0,0},
  {0,0,0,0,0,1,0,0,0,0},
  {0,0,0,0,0,0,1,0,0,0},
  {0,0,0,0,0,0,0,1,0,0},
  {0,0,0,0,0,0,0,0,1,0},
  {0,0,0,0,0,0,0,0,0,1}};
int active_dim[3] = {0,1,2};
int old_x,old_y;

void draw(int color)
{
  int x,y;
  int i,j,k;
  double v[n],w[n];


  for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++) v[j] = (data[INDEX(i,j)] - (max[j] + min[j]) / 2)/(max[j] - min[j]);
      for(j=0;j<n;j++) 
	{
	  w[j] = 0;
	  for(k=0;k<10;k++) w[j] += A[j][k] * v[k];
	}
      
      x = w[0] * X_SCALE + SCREEN_WIDTH / 2;
      y = w[1] * Y_SCALE + SCREEN_HEIGHT / 2;
      if (x >= 0 && y >= 0 && x < SCREEN_WIDTH && y < SCREEN_HEIGHT) point(x,y) = (color) ? pal[col[i]] : 0;
    }
}

void rotate(double theta, double phi)
{
  double b[2];
  int i;
  
  for(i=0;i<10;i++)
    {
      b[0] = A[i][active_dim[0]];
      b[1] = A[i][active_dim[1]];

      A[i][active_dim[0]] *= cos(theta);
      A[i][active_dim[1]] *= cos(theta);
      
      A[i][active_dim[0]] += sin(theta) * b[1];
      A[i][active_dim[1]] -= sin(theta) * b[0];
    }

  for(i=0;i<10;i++)
    {
      b[0] = A[i][active_dim[1]];
      b[1] = A[i][active_dim[2]];

      A[i][active_dim[1]] *= cos(phi);
      A[i][active_dim[2]] *= cos(phi);
      
      A[i][active_dim[1]] += sin(phi) * b[1];
      A[i][active_dim[2]] -= sin(phi) * b[0];
    }
  
}

int main(int argc,char ** argv)
{
  int flag;
  int i,j;
  FILE * fp;
  double theta,phi;

  SDL_Event event;

  screen_init();

  if (argc < 2) 
    {
      printf("Usage: plotn matrix_file [color_file]\n");
      exit(1);
    }
  
  if ((fp = fopen(argv[1],"r"))==NULL) {printf("Error opening file %s\n",argv[1]);exit(1);}
  n = cols(fp);

  if ((data = malloc(sizeof(double) * MAX_ROWS * n)) == NULL) {printf("Out of memory\n");exit(1);} 

  for(m=0;;m++)
    {
      if (m>=MAX_ROWS) {printf("Too many rows\n");exit(1);}

      for(i=0;i<n;i++)
	{
	  if (fscanf(fp,"%lf",&data[INDEX(m,i)])!=1) break;
	}
      if (i!=n) break;		    
    }
  
  if ((col = malloc(sizeof(unsigned) * m)) == NULL) {printf("Out of memory\n");exit(1);}

  for(j=0;j<n;j++) max[j] = min[j] = data[INDEX(0,j)];
  for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++)
	{
	  if (data[INDEX(i,j)] > max[j]) max[j] = data[INDEX(i,j)];
	  if (data[INDEX(i,j)] < min[j]) min[j] = data[INDEX(i,j)];
	}
    }
  fclose(fp);


  if (argc >= 3)
    {
      if ((fp = fopen(argv[2],"r"))==NULL) {printf("Error opening file %s\n",argv[2]);exit(1);}
      for(i=0;i<m;i++)
	{
	  if (fscanf(fp,"%d",&col[i])!=1) {printf("Error reading all colors\n");exit(1);}
	}
      fclose(fp);
    }
  else
    {
      for(i=0;i<m;i++) col[i] = 0;
    }

  flag = 1;
  draw(PEN_COLOR);
  refresh();



  old_x = -1;
  old_y = -1;
  while(flag)
    {
      //Redraw

      //Event loop
      if (SDL_WaitEvent(&event))
	{
	  switch(event.type)
	    {
	    case SDL_KEYDOWN:
	      if (event.key.keysym.sym == SDLK_q)
		{	         
		  flag = 0;
		}
	      else if (event.key.keysym.sym >= SDLK_0 && event.key.keysym.sym <= SDLK_9)
		{
		  i = event.key.keysym.sym - SDLK_0;
		  if (i != active_dim[1] && i!=active_dim[2])
		    {
		      active_dim[0] = active_dim[1];
		      active_dim[1] = active_dim[2];
		      active_dim[2] = i;
		    }
		}
	      else if (event.key.keysym.sym == SDLK_MINUS)
		{
		  draw(BACK_COLOR);
		  X_SCALE /= SCALE_RATIO;
		  Y_SCALE /= SCALE_RATIO;
		  draw(PEN_COLOR);
		  refresh();
		}
	      else if (event.key.keysym.sym == SDLK_EQUALS)
		{
		  draw(BACK_COLOR);
		  X_SCALE *= SCALE_RATIO;
		  Y_SCALE *= SCALE_RATIO;
		  draw(PEN_COLOR);
		  refresh();
		}
	      break;
	    case SDL_MOUSEMOTION:
	      while(SDL_PollEvent(&event));
	      if (old_x == -1)
		{
		  old_x = event.motion.x;
		  old_y = event.motion.y;
		  break;
		}	     
	      draw(BACK_COLOR);
	      theta = (event.motion.x - old_x) * TWO_PI / SCREEN_WIDTH;
	      phi   = (event.motion.y - old_y) * TWO_PI / SCREEN_HEIGHT;
	      rotate(theta,phi);
	      draw(PEN_COLOR);
	      refresh();

	      old_x = event.motion.x;
	      old_y = event.motion.y;
	      break;
	      //	    case SDL_MOUSEBUTTONDOWN:
	      //	      switch(event.button.button)
	      //{
	      //case SDL_BUTTON_LEFT:
	      //  break;
	      //}
	      //	      break;
	    case SDL_QUIT:
	      flag = 0;
	      break;	
	    }	  
	}
    }

}  
