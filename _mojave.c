#include <SDL2/SDL.h>
#include <SDL2/SDL_blendmode.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_keyboard.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_render.h>
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL_timer.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define FPS 30
#define FRAME_DELAY (1000 / FPS)

#define IMAGE_BRUSH_FILE "/home/kevin/big-yellow-duck/mojave_brush.bmp"
#define IMAGE_PALETTE_FILE "/home/kevin/big-yellow-duck/mojave_palette.bmp"
#define IMAGE_ROTATION_FILE "/home/kevin/big-yellow-duck/mojave_rotation.bmp"

#define POINT_SCREEN 0
#define CONTROL_SCREEN 1
#define BRUSH_SCREEN 2
#define SCREENS 3
#define SCREEN_WIDTHS {800,400,800}
#define SCREEN_HEIGHTS {800,1000,200}
#define SCREEN_XPOSES {0,850,0}
#define SCREEN_YPOSES {0,0,900}

#define CONTROL_BG_COLOR 0x404040
#define CONTROL_FG_COLOR 0x3080ff
#define CONTROL_RADIUS 32
#define CONTROL_Y_STEP 80
#define CONTROL_LINE_STEP 0.01
#define CONTROL_ARROW_COLOR 0xffffff
#define CONTROL_ARROW_THICKNESS 1
#define CONTROL_BOX_RADIUS 32
#define CONTROL_NUM_BOX 3
#define CONTROL_BOX_FG 0xff0000
#define CONTROL_BOX_BG 0x0080ff
#define CONTROL_BOX_SELECT 0x000050
#define CONTROL_BOX_SHIFT 1.2
#define POINT_ZOOM 300
#define POINT_ZOOM_MULT 1.05
#define ROTATION_SPEED_MULT 1.05
#define ZOOM_GROVE_Y 100
#define ZOOM_GROVE_X 380
#define ZOOM_GROVE_WIDTH 5
#define ZOOM_GROVE_HEIGHT 400
#define ZOOM_SLIDER_WIDTH 30
#define ZOOM_SLIDER_HEIGHT 10
#define ZOOM_GROVE_COLOR 0x000000
#define ZOOM_SLIDER_COLOR 0xffffff
#define ZOOM_LOG_RATIO 40.0
#define SPEED_GROVE_X 340
#define ROTATION_MODE_OFF 0
#define ROTATION_MODE_MEDIUM 1
#define ROTATION_MODE_HIGH 2
#define ROTATION_MODE_NUM 3
#define RANDOM_SEED (lrand48())

#define BRUSH_BG 0x404040
#define BIT_RADIUS 8
#define BIT_X_SKIP 25
#define BIT_COLOR 0xffffff
#define BIT_SELECT_COLOR 0x8080ff
#define PALETTE_BOX_SIZE 30
#define PALETTE_BOX_SKIP 60
#define SELECTED_PALETTE_BOX_BG 0x000000
#define BRUSH_MODE_MARGIN 30
#define BRUSH_MODE_BOX_SIZE 50
#define ROTATION_MODE_MARGIN 10
#define ROTATION_MODE_BOX_SIZE 50
#define PALETTE_BOX_MARGIN 10
#define BIT_USED_COLOR 0xffff00
#define BIT_UNUSED_COLOR 0x909090
#define BRUSH_COLOR_MODES 2
#define BRUSH_COLOR_MODE_DIRECT 0
#define BRUSH_COLOR_MODE_HASH 1
#define BRUSH_SIZE_RATIO 1.05
#define DEFAULT_BRUSH_XSIZE 32
#define DEFAULT_BRUSH_YSIZE 32
#define OFFSCREEN -100

#define R_POWERS 13
#define RX_THETA_MAX 0.01
#define RY_THETA_MAX 0.01
#define KEYBOARD_ROTATION_DX 0
#define KEYBOARD_ROTATION_DY 10
#define UNDO_SIZE 1024

#define SQR(x) ((x)*(x))
#define LCG(x) ((134775813 * (x) + 2531011) & 0xffffff)
#define COLOR_HASH(x) LCG(LCG(x))
#define RANDOM(n) ((int)(drand48() * (n)))

#define point(i, x, y) pnt[i][(int)(x) + (int)(y)*SCREEN_WIDTH[i]]
#define AA(i, j) ((i)*dim + (j))
#define RR(k, i, j) ((k) * dim * dim + (i) * dim + (j))

// SDL window stuff
SDL_Window *screen[SCREENS];
SDL_Renderer * renderer[SCREENS];
SDL_Texture * texture[SCREENS];
unsigned * pnt[SCREENS];
int SCREEN_WIDTH[SCREENS] = SCREEN_WIDTHS;
int SCREEN_HEIGHT[SCREENS] = SCREEN_HEIGHTS;
int SCREEN_XPOS[SCREENS] = SCREEN_XPOSES;
int SCREEN_YPOS[SCREENS] = SCREEN_YPOSES;

// Global transform data
// A shape (dim, dim) is the projection onto the first 2 dimensions.
// R is such that R[k] = R^(2^k)
double *A;
double rotation_direction_exists = 0;
double *Rx;
double *Rx_inv;
double *Ry;
double *Ry_inv;
int dim;

// Global image stuff
SDL_Surface *image_brush;
SDL_Surface *image_palette;
SDL_Surface *image_rotation;

// box is shape (dim, num-boxes)
int (*box)[CONTROL_NUM_BOX];

// Pallete = (x >> mask_location) & 7
int mask_location = 0;

// Brush color
int selected_color = 0;

// Pallete mode modifier.
int use_color_hash = 0;

// brush info
int brush_color_mode = 0;
int brush_xsize = DEFAULT_BRUSH_XSIZE;
int brush_ysize = DEFAULT_BRUSH_YSIZE;
int brush_x = OFFSCREEN;
int brush_y = OFFSCREEN;
unsigned brush_color[8] = {
			   0xffffff,
			   0x0000ff,
			   0x00ff00,
			   0x00ffff,
			   0xff0000,
			   0xff00ff,
			   0xffff00,
			   0x808080};

// Rotation info
unsigned rotation_seed = 0;


// Sliders
double zoom_ratio = 1.0;
double rotation_speed = 1.0;

// Mode flags
int rotation_mode = ROTATION_MODE_OFF;
int rotation_mode_color[3] = {0xff0000, 0xffff00, 0x00ff00}; 
int brush_mode_on = 0;

// SDL refresh
void refresh(int i)
{
  SDL_UpdateTexture(texture[i], NULL, pnt[i], SCREEN_WIDTH[i] * sizeof(unsigned));
  SDL_RenderClear(renderer[i]);
  SDL_RenderCopy(renderer[i], texture[i], NULL, NULL);
  SDL_RenderPresent(renderer[i]);
}

// SDL screen init
void screen_init(int i, int init, char * name, int xpos, int ypos)
{
  if (init)
    {
      SDL_Init(SDL_INIT_VIDEO);
      atexit(SDL_Quit);
      
      screen[i] = SDL_CreateWindow(name, xpos, ypos,
				   SCREEN_WIDTH[i], SCREEN_HEIGHT[i],
				   SDL_WINDOW_RESIZABLE);
      renderer[i] = SDL_CreateRenderer(screen[i], -1, 0);
    }
  SDL_RenderClear(renderer[i]);
  SDL_GL_SwapWindow(screen[i]);    
  texture[i] = SDL_CreateTexture(renderer[i],SDL_PIXELFORMAT_ARGB8888,
			      SDL_TEXTUREACCESS_STREAMING,
			      SCREEN_WIDTH[i], SCREEN_HEIGHT[i]);
  SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY,
	      "linear");  // make the scaled rendering look smoother.
  SDL_RenderSetLogicalSize(renderer[i], SCREEN_WIDTH[i], SCREEN_HEIGHT[i]);

  if (pnt[i]) free(pnt[i]);
  if ((pnt[i] = malloc(SCREEN_WIDTH[i] * SCREEN_HEIGHT[i] * sizeof(unsigned))) == 0)
    {fprintf(stderr, "OUT OF MEMORY");exit(1);}
}

void blt(SDL_Surface * surface, int screen_index, int x, int y,
	 int width, int height, unsigned icon_color)
{
  for(int dy = 0; dy < height; dy++)
    for(int dx = 0; dx < width; dx++)
      {
	unsigned rx = surface->w * dx / width;
       	unsigned ry = surface->h * dy / height;
	unsigned c = * (uint32_t *) ((uint8_t *) surface->pixels
				     + ry * surface->pitch
				     + rx * surface->format->BytesPerPixel);
	if ((c & 0xffffff) == 0xffffff) continue;
	point(screen_index, x + dx, y + dy) = icon_color;
      }
}

// Get RGB color value from logical color label
unsigned get_color(unsigned color_value)
{
  if (brush_color_mode == BRUSH_COLOR_MODE_DIRECT)
    return brush_color[(color_value >> mask_location) & 7];
  else
    return COLOR_HASH(color_value);
}

// The global transformation to point window coordinates (define by A)
void transform(double * data_point, double * out_x, double * out_y)
{
  double x = 0.0;
  double y = 0.0;
  for(int j=0; j < dim; j++)
    {
      x += A[AA(j,0)] * data_point[j];
      y += A[AA(j,1)] * data_point[j];
    }
  double x0 = SCREEN_WIDTH[POINT_SCREEN]/2.0;
  double y0 = SCREEN_HEIGHT[POINT_SCREEN]/2.0;
  double xs = x * POINT_ZOOM * zoom_ratio + x0;
  double ys = y * POINT_ZOOM * zoom_ratio + y0;

  // Use clipped points
  if (xs < 0) xs = 0.0;
  if (ys < 0) ys = 0.0;
  if (xs >= SCREEN_WIDTH[POINT_SCREEN]) xs = SCREEN_WIDTH[POINT_SCREEN] - 1;
  if (ys >= SCREEN_HEIGHT[POINT_SCREEN]) ys = SCREEN_HEIGHT[POINT_SCREEN] - 1;

  *out_x = xs;
  *out_y = ys;
}

// Draws one of the sliders onto the control window.
void draw_slider(double log_ratio, double ratio, double grove_width,
		 double grove_height, double grove_x, double grove_y,
		 double slider_width, double slider_height,
		 unsigned grove_color, unsigned slider_color)
{
  int slider_y = (log_ratio * log2(ratio)) + grove_height / 2.0;
  if (slider_y < 0) slider_y = 0;
  if (slider_y >= grove_height) slider_y = grove_height - 1;
  for(int y = grove_y; y < grove_y + grove_height; y++)
    for(int x = grove_x; x < grove_x + grove_width; x++)
      point(CONTROL_SCREEN, x, y) = grove_color;
  for(int dy = -slider_height/2; dy <= slider_height/2; dy++)
    for(int dx = -slider_width/2; dx <= slider_width/2; dx++)
      point(CONTROL_SCREEN, grove_x + grove_width/2 + dx,
	    grove_y + dy + slider_y) = slider_color;
}

// Draws the control window
void draw_controls()
{
  // Background fill
  for(int y=0; y<SCREEN_HEIGHT[1]; y++)
    for(int x=0; x<SCREEN_WIDTH[1]; x++)
      point(CONTROL_SCREEN,x,y) = CONTROL_BG_COLOR;

  // Main controls
  for(int i = 0; i < dim; i++)
    {
      // Circle displays
      for(int y=0;y<2*CONTROL_RADIUS;y++)
	for(int x=0;x<2*CONTROL_RADIUS;x++)
	  {
	    int x0 = x - CONTROL_RADIUS;
	    int y0 = y - CONTROL_RADIUS;
	    if (SQR(x0) + SQR(y0) <= SQR(CONTROL_RADIUS))
	      point(CONTROL_SCREEN,x,y+i*CONTROL_Y_STEP) = CONTROL_FG_COLOR;
	  }     
      double dx = CONTROL_RADIUS * A[AA(i,0)];
      double dy = CONTROL_RADIUS * A[AA(i,1)];
      for(double r=0.0; r<1.0; r+=CONTROL_LINE_STEP)
	{
	  for(int dx2=-CONTROL_ARROW_THICKNESS;dx2<=CONTROL_ARROW_THICKNESS;dx2++)
	    for(int dy2=-CONTROL_ARROW_THICKNESS;dy2<=CONTROL_ARROW_THICKNESS;dy2++)
	      point(CONTROL_SCREEN,
		    dx2 + (int)(CONTROL_RADIUS+dx*r),
		    dy2 + (int)(i*CONTROL_Y_STEP+CONTROL_RADIUS+dy*r)
		    ) = CONTROL_ARROW_COLOR;
	}

      // Boxes
      for(int bi = 0; bi < CONTROL_NUM_BOX; bi++)
	for(int y=-CONTROL_BOX_RADIUS;y<=CONTROL_BOX_RADIUS;y++)
	  for(int x=-CONTROL_BOX_RADIUS;x<=CONTROL_BOX_RADIUS;x++)
	    {
	      int color;
	      if (y==-CONTROL_BOX_RADIUS || x==-CONTROL_BOX_RADIUS ||\
		  y==CONTROL_BOX_RADIUS || x==CONTROL_BOX_RADIUS)
		color = CONTROL_BOX_FG;
	      else
		color = box[i][bi] ? CONTROL_BOX_SELECT : CONTROL_BOX_BG;
	      point(CONTROL_SCREEN,
		    x+CONTROL_Y_STEP*(bi+CONTROL_BOX_SHIFT) + CONTROL_RADIUS,
		    y+i*CONTROL_Y_STEP + CONTROL_RADIUS
		    ) = color;
	    }      
    }

  // Rotation mode indicator
  blt(image_rotation, CONTROL_SCREEN,
      SCREEN_WIDTH[CONTROL_SCREEN] - ROTATION_MODE_MARGIN - ROTATION_MODE_BOX_SIZE,
      ROTATION_MODE_MARGIN,
      ROTATION_MODE_BOX_SIZE, ROTATION_MODE_BOX_SIZE,
      rotation_mode_color[rotation_mode]);

  // Draw sliders
  draw_slider(ZOOM_LOG_RATIO, zoom_ratio, ZOOM_GROVE_WIDTH, ZOOM_GROVE_HEIGHT,
	      ZOOM_GROVE_X, ZOOM_GROVE_Y, ZOOM_SLIDER_WIDTH, ZOOM_SLIDER_HEIGHT,
	      ZOOM_GROVE_COLOR, ZOOM_SLIDER_COLOR);
  draw_slider(ZOOM_LOG_RATIO, rotation_speed, ZOOM_GROVE_WIDTH, ZOOM_GROVE_HEIGHT,
	      SPEED_GROVE_X, ZOOM_GROVE_Y, ZOOM_SLIDER_WIDTH, ZOOM_SLIDER_HEIGHT,
	      ZOOM_GROVE_COLOR, ZOOM_SLIDER_COLOR);
}

// Draws the brush window
void draw_palette(int num_data, double (*data)[dim], int32_t * color)
{
  unsigned mask = 0;
  for(int i = 0; i < num_data; i++) mask |= color[i];
  
  for(int y=0;y<SCREEN_HEIGHT[BRUSH_SCREEN];y++)
    for(int x=0;x<SCREEN_WIDTH[BRUSH_SCREEN];x++)
      point(BRUSH_SCREEN,x,y) = BRUSH_BG;

  for(int bx=0;bx<30;bx++)
    {
      int bit_color = ((mask >> bx) & 1) ? BIT_USED_COLOR : BIT_UNUSED_COLOR;
      for(int x=-BIT_RADIUS;x<=BIT_RADIUS;x++)
	for(int y=-BIT_RADIUS;y<=BIT_RADIUS;y++)
	  {
	    if (SQR(x) + SQR(y) <= SQR(BIT_RADIUS))
	      {
		point(BRUSH_SCREEN, x + (bx + 1) * BIT_X_SKIP,
		      30 + y + BIT_RADIUS) = bit_color;
	      }
	  }

      if (bx % 3 != 1) continue;

      bit_color = ((mask >> bx) & 7) ? BIT_USED_COLOR : BIT_UNUSED_COLOR;

      for(int x=-BIT_RADIUS;x<=BIT_RADIUS;x++)
	for(int y=-BIT_RADIUS;y<=BIT_RADIUS;y++)
	  {
	    if (SQR(x) + SQR(y) <= SQR(BIT_RADIUS))
	      {
		point(BRUSH_SCREEN, x + (bx + 1) * BIT_X_SKIP,
		      30 + y + BIT_RADIUS + BIT_X_SKIP + 10) = bit_color;
	      }
	  }
    }

  for(int x = 0; x < 3 * BIT_X_SKIP; x++)
    for(int y = 0; y < 5; y++)
      {
	point(BRUSH_SCREEN, BIT_X_SKIP*mask_location + x+10,  y + BIT_X_SKIP / 2)
	  = BIT_SELECT_COLOR;
	point(BRUSH_SCREEN, BIT_X_SKIP*mask_location + x+10,  y + 30 + BIT_X_SKIP)
	  = BIT_SELECT_COLOR;
      }

  for(int c=0; c < 8; c++)
    {
      if (c == selected_color)
	{
	  for(int x=-PALETTE_BOX_MARGIN;x<PALETTE_BOX_SIZE+PALETTE_BOX_MARGIN;x++)
	    for(int y=-PALETTE_BOX_MARGIN;y<PALETTE_BOX_SIZE+PALETTE_BOX_MARGIN;y++)
	      point(BRUSH_SCREEN, BIT_X_SKIP/2 + x + c*PALETTE_BOX_SKIP,
		    y+130) = SELECTED_PALETTE_BOX_BG;      
	}
      
      for(int x=0;x<PALETTE_BOX_SIZE;x++)
	for(int y=0;y<PALETTE_BOX_SIZE;y++)
	  point(BRUSH_SCREEN, BIT_X_SKIP/2 + x + c*PALETTE_BOX_SKIP,
		y+130) = brush_color[c];      
    }

  // Draw brush mode icon
  blt(image_brush, BRUSH_SCREEN,
      SCREEN_WIDTH[BRUSH_SCREEN] - BRUSH_MODE_MARGIN - BRUSH_MODE_BOX_SIZE,
      SCREEN_HEIGHT[BRUSH_SCREEN] - BRUSH_MODE_MARGIN - BRUSH_MODE_BOX_SIZE,
      BRUSH_MODE_BOX_SIZE, BRUSH_MODE_BOX_SIZE,
      brush_mode_on ? 0x00ff00 : 0xff0000);

  // Draw color mode icon
  blt(image_palette, BRUSH_SCREEN,
      SCREEN_WIDTH[BRUSH_SCREEN] - 2 * BRUSH_MODE_MARGIN - 2 * BRUSH_MODE_BOX_SIZE,
      SCREEN_HEIGHT[BRUSH_SCREEN] - BRUSH_MODE_MARGIN - BRUSH_MODE_BOX_SIZE,
      BRUSH_MODE_BOX_SIZE, BRUSH_MODE_BOX_SIZE,
      brush_color_mode ? 0x00ff00 : 0xff0000);
}

// Draws the main view - points window
void draw_points(int num_data, double (*data)[dim], int32_t * color)
{  
  // Draw points
  for(int i=0; i < num_data; i++)
    {
      double x,y;
      transform(data[i], &x, &y);
      point(POINT_SCREEN,x,y) = get_color(color[i]);
    }

  if (brush_xsize == 0 || brush_ysize == 0) return;
  
  // Draw horizontal lines of brush
  int xstep = (brush_xsize >= 0) ? 1 : -1;
  for(int dx = 0; dx != brush_xsize; dx+=xstep)
    {
      int x = brush_x + dx;
      int y = brush_y;
      if (x>=0 && y>=0 && x <= SCREEN_WIDTH[POINT_SCREEN] &&
	  y <= SCREEN_HEIGHT[POINT_SCREEN])
	point(POINT_SCREEN, x, y) = brush_color[selected_color];
      y = brush_y + brush_ysize;
      if (x>=0 && y>=0 && x <= SCREEN_WIDTH[POINT_SCREEN] &&
	  y <= SCREEN_HEIGHT[POINT_SCREEN])
	point(POINT_SCREEN, x, y) = brush_color[selected_color];
    }

  // Draw vertical lines of brush
  int ystep = (brush_ysize >= 0) ? 1 : -1;
  for(int dy = 0; dy != brush_ysize; dy+=ystep)
    {
      int y = brush_y + dy;
      int x = brush_x;
      if (x>=0 && y>=0 && x <= SCREEN_WIDTH[POINT_SCREEN] &&
	  y <= SCREEN_HEIGHT[POINT_SCREEN])
	point(POINT_SCREEN, x, y) = brush_color[selected_color];
      x = brush_x + brush_xsize;
      if (x>=0 && y>=0 && x <= SCREEN_WIDTH[POINT_SCREEN] &&
	  y <= SCREEN_HEIGHT[POINT_SCREEN])
	point(POINT_SCREEN, x, y) = brush_color[selected_color];
    }

}

// Sets out to the identity
void SO_clear(double * out)
{
  for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
      out[AA(i,j)] = i==j;
}

// Multiplies rotations matrices
void SO_mult(double * out, double * A1, double * A2)
{
  double temp[SQR(dim)];
  for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
      {
	temp[AA(i,j)] = 0.0;
	for(int k=0;k<dim;k++) temp[AA(i, j)] += A1[AA(i,k)] * A2[AA(k,j)];
    }
  for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
      out[AA(i,j)] = temp[AA(i,j)];
}

// Inverts a rotation matrix (easy since it is just atranspose.)
void SO_inverse(double *out, double *in)
{
  for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
      out[AA(i,j)] = in[AA(j,i)];
}

// Create a rotation matrix whiche rotates in the dim1/dim2 plane
void SO_pair(double * out, int dim1, int dim2, double theta)
{
  for(int i=0;i<dim;i++)
    for(int j=0;j<dim;j++)
      out[AA(i,j)] = (i==j);
  out[AA(dim1,dim1)] = cos(theta);
  out[AA(dim1,dim2)] = -sin(theta);
  out[AA(dim2,dim1)] = sin(theta);
  out[AA(dim2,dim2)] = cos(theta);
}

// Create powers of a rotation matrix to expedite special powers.
void SO_powers(double * R)
{
  for(int i=1;i<R_POWERS;i++)
    SO_mult(&R[i*dim*dim], &R[0], &R[(i-1)*dim*dim]);
}

// Special power (using SO_powers)
void SO_power(double *out, double *R, int power)
{
  SO_clear(out);
  for(int i=0;i<R_POWERS;i++)
    {
      if ((power >> i) & 1) SO_mult(out,out,&R[i*dim*dim]);
    }
}

// Compute the dz-th power of Rz
void rotate_dim(double * Rz, double * Rz_inv, int dz)
{
  double Z[SQR(dim)];
  if (dz > 0)
    {
      SO_power(Z, Rz, dz);
      SO_mult(A, Z, A);
    }
  else if (dz < 0)
    {
      SO_power(Z, Rz_inv, -dz);
      SO_mult(A, Z, A);
    }
}

// Rotate using dx and dy (powers)
void SO_rotate(int dx, int dy)
{
  rotate_dim(Rx, Rx_inv, dx);
  rotate_dim(Ry, Ry_inv, dy);
}

// Pick new rotation directions for dx and dy
void new_rotation_direction(unsigned seed)
{
  rotation_seed = seed;
  srand48(seed);
  
  // Reset rotation directions
  int cnt = 0;
  for(int i=0;i<dim;i++) cnt += box[i][2];
  if (cnt < 2) return;     
  
  // Set up rotations Rx a pair
  int ind_1 = RANDOM(cnt);
  int ind_2 = RANDOM(cnt-1);
  int dim1, dim2;
  if (ind_2 >= ind_1) ind_2++;
  cnt = 0;
  for(int i=0;i<dim;i++)
    {
      if (cnt == ind_1) dim1 = i;
      if (cnt == ind_2) dim2 = i;
      cnt += box[i][2];
    }
  double theta = (1.0-2.0*drand48()) * RX_THETA_MAX * rotation_speed;
  SO_pair(&Rx[0], dim1, dim2, theta);
  SO_powers(Rx);
  SO_inverse(&Rx_inv[0], &Rx[0]); 
  SO_powers(Rx_inv);
  
  // Set up rotations Ry a bunch of pairs
  SO_clear(&Ry[0]);
  for(int k=0;k<SQR(dim);k++)
    {
      int ind_1 = RANDOM(cnt);
      int ind_2 = RANDOM(cnt-1);
      int dim1, dim2;
      if (ind_2 >= ind_1) ind_2++;
      cnt = 0;
      for(int i=0;i<dim;i++)
	{
	  if (cnt == ind_1) dim1 = i;
	  if (cnt == ind_2) dim2 = i;
	  cnt += box[i][2];
	}
      double theta = (1.0-2.0*drand48()) * RY_THETA_MAX * rotation_speed;
      double P[SQR(dim)];
      SO_pair(P, dim1, dim2, theta);
      SO_mult(&Ry[0], P, &Ry[0]);
    }      
  SO_powers(Ry);
  SO_inverse(&Ry_inv[0], &Ry[0]); 
  SO_powers(Ry_inv);
  
  // Set rotation direction flag
  rotation_direction_exists = 1;
}

void change_rotation_mode()
{
  brush_mode_on = 0;
  if (++rotation_mode >= ROTATION_MODE_NUM)
    rotation_mode = 0;
  new_rotation_direction(RANDOM_SEED);
}

// Move the slides given an (x,y) choice
void move_sliders(int x, int y)
{
  int middle_grove_x = (ZOOM_GROVE_X + SPEED_GROVE_X) / 2;
  double ratio = pow(2, (y - ZOOM_GROVE_Y - ZOOM_GROVE_HEIGHT / 2.0)
		     / ZOOM_LOG_RATIO);
  if (x < middle_grove_x)
    rotation_speed = ratio;
  else
    zoom_ratio = ratio;
}

// data must be normalized to be in [-1,1]
// color will be modified in place
void mojave(double * data_flat, int32_t * color, int num_data, int dim_in,
	    char * name)
{
  SDL_SetHint(SDL_HINT_MOUSE_FOCUS_CLICKTHROUGH, "1");
  
  dim = dim_in;  
  if (dim <= 1) return;
  double (*data)[dim] = (double (*)[dim]) data_flat;

  int32_t undo[UNDO_SIZE][num_data];
  int undo_length = 1;
  int max_undo_length = 1;
  for(int i=0;i<num_data;i++) undo[0][i] = color[i];
  
  // Set up initial transform (A) \in SO(dim)
  if ((A = malloc(SQR(dim) * sizeof(double))) == NULL)
    {fprintf(stderr, "Out of memory\n");exit(1);}  
  for(int y=0;y<dim;y++)
    for(int i=0;i<dim;i++)
      A[AA(y,i)] = (i==y);

  // Set up initial rotations (R,R_inv)
  if ((Rx = malloc(R_POWERS * SQR(dim) * sizeof(double))) == NULL)
    {fprintf(stderr, "Out of memory\n");exit(1);}  
  if ((Rx_inv = malloc(R_POWERS * SQR(dim) * sizeof(double))) == NULL)
    {fprintf(stderr, "Out of memory\n");exit(1);}
  if ((Ry = malloc(R_POWERS * SQR(dim) * sizeof(double))) == NULL)
    {fprintf(stderr, "Out of memory\n");exit(1);}  
  if ((Ry_inv = malloc(R_POWERS * SQR(dim) * sizeof(double))) == NULL)
    {fprintf(stderr, "Out of memory\n");exit(1);}
  for(int k=0;k<R_POWERS;k++)
    for(int y=0;y<dim;y++)
      for(int i=0;i<dim;i++)
	{
	  Rx_inv[RR(k,y,i)] = Rx[RR(k,y,i)] = (i==y);
	  Ry_inv[RR(k,y,i)] = Ry[RR(k,y,i)] = (i==y);
	}
  
  // Set up initial controls (boxes)
  if ((box = malloc(dim * CONTROL_NUM_BOX * sizeof(int))) == NULL)
    {fprintf(stderr, "Out of memory\n");exit(1);}
  for(int y=0; y<dim; y++)
    for(int i=0; i<CONTROL_NUM_BOX; i++)
      box[y][i] = 0;
  box[0][0] = 1;
  box[1][1] = 1;

  screen_init(POINT_SCREEN, 1, name, SCREEN_XPOS[POINT_SCREEN],
	      SCREEN_YPOS[POINT_SCREEN]);
  screen_init(CONTROL_SCREEN, 1, "Controls", SCREEN_XPOS[CONTROL_SCREEN],
	      SCREEN_YPOS[CONTROL_SCREEN]);
  screen_init(BRUSH_SCREEN, 1, "Brush", SCREEN_XPOS[BRUSH_SCREEN],
	      SCREEN_YPOS[BRUSH_SCREEN]);
  SDL_SetRenderDrawColor(renderer[POINT_SCREEN],0,0,0,255);

  image_brush = SDL_LoadBMP(IMAGE_BRUSH_FILE);
  image_palette = SDL_LoadBMP(IMAGE_PALETTE_FILE);
  image_rotation = SDL_LoadBMP(IMAGE_ROTATION_FILE);
  
  SDL_Event event, dummy;
  int flag = 1;
  int refresh_flag = 1;
  int last_mouse_x = -1;
  int last_mouse_y = -1;
  int mouse_x, mouse_y;
  while(flag)
    {
      unsigned frame_start = SDL_GetTicks();

      // Non-event driven rotation
      if (rotation_mode == ROTATION_MODE_HIGH)
	{
	  SO_rotate(KEYBOARD_ROTATION_DX, KEYBOARD_ROTATION_DY);
	  refresh_flag = 1;
	}
      
      // Refresh logic
      if (refresh_flag)
	{
	  // Point Screen
	  for(int x=0;x<SCREEN_WIDTH[POINT_SCREEN];x++)
	    for(int y=0;y<SCREEN_WIDTH[POINT_SCREEN];y++) point(POINT_SCREEN,x,y)=0;
	  draw_points(num_data, data, color);
	  refresh(POINT_SCREEN);

	  // Control screen
	  draw_controls();
	  refresh(CONTROL_SCREEN);

	  // Brush screen
	  draw_palette(num_data, data, color);
	  refresh(BRUSH_SCREEN);
	}
      refresh_flag = 0;

      // Event loop
      if (SDL_PollEvent(&event))
	{
	  while (SDL_PollEvent(&dummy));

	  switch(event.type)
	    {
            case SDL_WINDOWEVENT:	      
	      if (event.window.event == SDL_WINDOWEVENT_RESIZED)
		{
		  if (event.window.windowID ==
		      SDL_GetWindowID(screen[POINT_SCREEN]))
		    {		      
		      SCREEN_WIDTH[POINT_SCREEN] = event.window.data1;
		      SCREEN_HEIGHT[POINT_SCREEN] = event.window.data2;
		      screen_init(POINT_SCREEN, 0, name, 0, 0);
		      refresh_flag = 1;
		    }
		}
	      if (event.window.event == SDL_WINDOWEVENT_CLOSE)
		{
		  flag = 0;
		}
	      break;
	    case SDL_KEYDOWN:
	      switch(event.key.keysym.sym)
		{
		case SDLK_q:
		  flag = 0;
		  break;
		case SDLK_c:
		  if (++brush_color_mode >= BRUSH_COLOR_MODES)  brush_color_mode=0;
		  refresh_flag = 1;
		  break;
		case SDLK_SPACE:
		  brush_x = OFFSCREEN;
		  brush_y = OFFSCREEN;
		  brush_xsize = DEFAULT_BRUSH_XSIZE;
		  brush_ysize = DEFAULT_BRUSH_YSIZE;
		  rotation_direction_exists = 0;
		  new_rotation_direction(RANDOM_SEED);
		  brush_mode_on = 0;
		  brush_color_mode = 0;
		  refresh_flag = 1;
		  break;
		case SDLK_s:
		  brush_x = OFFSCREEN;
		  brush_y = OFFSCREEN;
		  brush_xsize = DEFAULT_BRUSH_XSIZE;
		  brush_ysize = DEFAULT_BRUSH_YSIZE;
		  SO_clear(A);
		  rotation_direction_exists = 0;
		  rotation_mode = ROTATION_MODE_OFF;
		  brush_mode_on = 0;
		  brush_color_mode = 0;
		  zoom_ratio = 1.0;
		  refresh_flag = 1;
		  break;
		case SDLK_r:
		  change_rotation_mode();
		  refresh_flag = 1;
		  break;
		case SDLK_b:
		  rotation_mode = ROTATION_MODE_OFF;
		  brush_mode_on = !brush_mode_on;
		  refresh_flag = 1;
		  break;
		case SDLK_o:
		  SDL_GetMouseState(&mouse_x, &mouse_y);
		  if (event.window.windowID ==
		      SDL_GetWindowID(screen[POINT_SCREEN]))
		    {
		      double min_dist_sqr = INFINITY;
		      unsigned best_color = 0x000000;
		      for(int i = 0; i < num_data; i++)
			{
			  double trans_x, trans_y;
			  transform(data[i], &trans_x, &trans_y);
			  double this_dist_sqr = SQR(trans_x - mouse_x) +
			    SQR(trans_y - mouse_y);
			  if (this_dist_sqr < min_dist_sqr)
			    {
			      min_dist_sqr = this_dist_sqr;
			      best_color = color[i];
			    }
			}
		      selected_color = 0; //best_color; @@@@@@@@@@@@
		    }
		  break;
		case SDLK_z:
		  if (undo_length > 1)
		    {
		      undo_length--;
		      for(int i=0;i<num_data;i++) color[i] = undo[undo_length-1][i];
		      refresh_flag = 1;
		    }
		  break;
		case SDLK_y:
		  if (undo_length < UNDO_SIZE && undo_length < max_undo_length)
		    {
		      for(int i=0;i<num_data;i++) color[i] = undo[undo_length][i];
		      undo_length++;
		      refresh_flag = 1;
		    }
		  break;
		case SDLK_PERIOD:
		  refresh_flag = 1;
                  SO_rotate(KEYBOARD_ROTATION_DX,
			    KEYBOARD_ROTATION_DY);
                  break;
		case SDLK_COMMA:
		  refresh_flag = 1;
		  SO_rotate(-KEYBOARD_ROTATION_DX,
			    -KEYBOARD_ROTATION_DY);
		  break;
		case SDLK_MINUS:
		  zoom_ratio /= POINT_ZOOM_MULT;
		  refresh_flag = 1;
		  break;
		case SDLK_EQUALS:
		  zoom_ratio *= POINT_ZOOM_MULT;
		  refresh_flag = 1;
		  break;
		case SDLK_LEFTBRACKET:
		  rotation_speed /= ROTATION_SPEED_MULT;
		  new_rotation_direction(rotation_seed);
		  refresh_flag = 1;
		  break;
		case SDLK_RIGHTBRACKET:
		  rotation_speed *= ROTATION_SPEED_MULT;
		  new_rotation_direction(rotation_seed);
		  refresh_flag = 1;
		  break;		 
		}
	    case SDL_MOUSEMOTION:
	      SDL_GetMouseState(&mouse_x, &mouse_y);
 	      if (event.window.windowID ==
		  SDL_GetWindowID(screen[POINT_SCREEN]))
		{
		  if (brush_mode_on || (event.motion.state & SDL_BUTTON_LMASK))
		    {
		      if (SDL_GetModState() & KMOD_CTRL)
			{
			  brush_xsize = mouse_x - brush_x;
			  brush_ysize = mouse_y - brush_y;
			  refresh_flag = 1;
			}
		      else
			{
			  brush_x = mouse_x - brush_xsize;
			  brush_y = mouse_y - brush_ysize;
			  for(int i=0;i<num_data;i++)
			    {
			      double x,y;
			      transform(data[i], &x, &y);
			      int x1 = (brush_xsize >= 0)
				? brush_x : brush_x + brush_xsize;
			      int x2 = (brush_xsize < 0)
				? brush_x : brush_x + brush_xsize;
			      int y1 = (brush_ysize >= 0)
				? brush_y : brush_y + brush_ysize;
			      int y2 = (brush_ysize < 0)
				? brush_y : brush_y + brush_ysize;
			      if (x >= x1 && x < x2 && y >= y1 && y < y2)
				color[i] = (color[i] &
					    ((7 << mask_location) ^ 0xffffffff)) ^
				  (selected_color << mask_location);
			    }
			  refresh_flag = 1;
			}
		    }
		  if (rotation_mode == ROTATION_MODE_MEDIUM
		    || (event.motion.state & SDL_BUTTON_RMASK))
		    {
		      brush_x = OFFSCREEN;
		      brush_y = OFFSCREEN;
		      
		      // Advance rotation
		      if (!rotation_direction_exists)
			new_rotation_direction(RANDOM_SEED);
		      if (last_mouse_x >= 0 && last_mouse_y >= 0)
			{
			  int dx = mouse_x - last_mouse_x;
			  int dy = mouse_y - last_mouse_y;
			  SO_rotate(dx,dy);
			}		      
		      last_mouse_x = mouse_x;
		      last_mouse_y = mouse_y;
		      refresh_flag = 1;
		    }
		  else
		    {
		      last_mouse_x = -1;
		      last_mouse_y = -1;
		    }
		}
	      else if (event.window.windowID ==
		       SDL_GetWindowID(screen[CONTROL_SCREEN]) &&
		       event.button.y >= ZOOM_GROVE_Y &&
		       event.button.y < ZOOM_GROVE_Y + ZOOM_GROVE_HEIGHT &&
		       event.button.button == SDL_BUTTON_LEFT)
		{
		  move_sliders(mouse_x, mouse_y);
		  refresh_flag = 1;
		}
	      break;	      
	    case SDL_MOUSEBUTTONDOWN:
	      switch(event.button.button)
		{
		case SDL_BUTTON_LEFT:		  
		  if (event.window.windowID ==
		      SDL_GetWindowID(screen[CONTROL_SCREEN]))
		    {
		      if (event.button.x < SPEED_GROVE_X - ZOOM_SLIDER_WIDTH/2.0)
			{			
			  int bi = round(((double) event.button.x - CONTROL_RADIUS)
					 / CONTROL_Y_STEP - CONTROL_BOX_SHIFT);
			  double xx = event.button.x -
			    CONTROL_Y_STEP*(bi+CONTROL_BOX_SHIFT) - CONTROL_RADIUS;
			  int i = round(((double) event.button.y - CONTROL_RADIUS)
					/ CONTROL_Y_STEP);
			  double yy = event.button.y - CONTROL_Y_STEP*i -
			    CONTROL_RADIUS;
			  if (bi>=0 && bi < CONTROL_NUM_BOX && i>=0 && i < dim &&
			      xx >= -CONTROL_BOX_RADIUS &&
			      xx <= CONTROL_BOX_RADIUS &&
			      yy >= -CONTROL_BOX_RADIUS && yy <= CONTROL_BOX_RADIUS)
			    {
			      if (bi >= 0 && bi < 2)
				{
				  int flag = 0;
				  for(int j=0;j<dim;j++)
				    {
				      flag |= box[j][2];
				      box[j][2] = 0;
				    }
				  if (flag) SO_clear(A);
				  if (!box[i][bi] && !box[i][!bi])
				    {
				      box[i][bi] = 1;
				      A[AA(i,bi)] = 1.0;
				      for(int j=0;j<dim;j++)
					{
					  if (i!=j)
					    {
					      A[AA(j,bi)] = 0;
					      box[j][bi] = 0;
					    }
					}
				    }				
				}
			      if (bi==2)
				{
				  int flag = 0;
				  for(int j=0;j<dim;j++)
				    {
				      flag |= box[j][0];
				      flag |= box[j][1];
				      box[j][0] = box[j][1] = 0;
				    }
				  if (flag) SO_clear(A);
				  rotation_direction_exists = 0;
				  box[i][bi] = !box[i][bi];
				}
			      refresh_flag = 1;
			    }
			  brush_x = OFFSCREEN;
			  brush_y = OFFSCREEN;
			}
		      else if (event.button.y >= ZOOM_GROVE_Y &&
			       event.button.y < ZOOM_GROVE_Y + ZOOM_GROVE_HEIGHT)
			{
			  move_sliders(event.button.x, event.button.y);
			  refresh_flag = 1;
			}
		      else if (event.button.y >= ROTATION_MODE_MARGIN &&
			       event.button.y < ROTATION_MODE_MARGIN +
			       ROTATION_MODE_BOX_SIZE &&
			       event.button.x >= SCREEN_WIDTH[CONTROL_SCREEN]
			       - ROTATION_MODE_MARGIN - ROTATION_MODE_BOX_SIZE &&
			       event.button.x <= SCREEN_WIDTH[CONTROL_SCREEN]
			       - ROTATION_MODE_MARGIN)
			{
			  change_rotation_mode();
			  refresh_flag = 1;
			}
		    }
		  if (event.window.windowID ==
		      SDL_GetWindowID(screen[BRUSH_SCREEN]))
		    {
		      if (event.button.y < 100)
			{
			  int bx = 3 * (event.button.x / (3 * BIT_X_SKIP));
			  if (bx < 0) bx = 0;
			  if (bx > 27) bx = 27;
			  mask_location = bx;			  
			  refresh_flag = 1;
			}
		      else
			{
			  int c = (event.button.x - BIT_X_SKIP/2)
			    / PALETTE_BOX_SKIP;
			  if (c < 0) c=0;
			  if (c > 7) c=7;
			  selected_color = c;
			  refresh_flag = 1;
			}
		    }
		  break;
		}
	      break;
	    case SDL_MOUSEBUTTONUP:
	      if (undo_length < UNDO_SIZE)
		{
		  int save_undo = 0;
		  if (undo_length == 0) save_undo = 1;
		  else
		    {
		      for(int i = 0; i < num_data; i++)
			{
			  if (color[i] != undo[undo_length - 1][i])
			    {
			      save_undo = 1;
			      break;
			    }
			}
		    }
		  if (save_undo)
		    {
		      for(int i=0;i<num_data;i++) undo[undo_length][i] = color[i];
		      undo_length++;
		      max_undo_length = undo_length;
		    }
		}
	      break;
	    }
	}

      
      unsigned frame_time = SDL_GetTicks() - frame_start;
      if (FRAME_DELAY > frame_time) SDL_Delay(FRAME_DELAY - frame_time);
    }

  SDL_Quit();
}

