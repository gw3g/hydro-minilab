/*
 * Author: greg jackson
 * Date: Dec 05 2015
 * (1+1)-dimensional advection
 * with: central diff
 *       upwind
 *       shasta
 */
#include "../include/shasta.h"

double x[N], w[N], u[N];
double x0;

double f0(double x) { 
  return 1./(1.+x*x); 
  /*return 1./(exp(x*2.)+1.);*/
}

/*
 * uncorrected transport: (u=const)
 */
void centre(double w[], double u) {

  for (int l=0; l<N; l++) {     w[l] -= u*dt*( w[ (l+1)%N ] - w[ (l-1+N)%N ] )/(2.*dx);      }

  return;
}

void upwind(double w[], double u) {

  for (int l=0; l<N; l++) {     w[l] -= u*dt*( w[ l ] - w[ (l-1+N)%N ] )/(dx);               }

  return;
}

/*
 * SHASTA [ Boris & Book, 1973 ]
 *
 */
int coord( int i, bnd X ) {                  // NB: only treats {-1, 0, ..., N-1, N}
  switch(X) {
    case P: return (i+N)%N;                  // periodic
    case F: return -(i)/N + (N-i-1)/N + i;   // fixed ends
  }
  return 0;                                  // checked: 16.01.06
}

void init() {
  /*
   *  Initial Conditions:   x    = { x_1, ... , x_N }   -  uniform grid
   *                        u    = { w_1, ... , u_N }   -  velocity field
   *                        w    = { w_1, ... , w_N }   -  scalar quantity
   */
  for (int i=0; i<N; i++) {
                            x[i] = x0 + dx*((double) i)   ; 
                            u[i] = 1.                   ;
                            /*w[i] = f0( x[i] - 5. )      ;*/
                            w[i] = f0( x[i] - 5. )      ;
  }
  return;
}

/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[40];

void eval_surf (int Nt, double del_t) { dt = del_t;

  init();
  sprintf(fname, "out/data/surf, (forward).dat", ( (double) Nt )*dt);
  file = fopen(fname, "w+");

  for (int n=0; n<Nt;n++) { 

    centre( w, 1.);       // central diff update
/* 
    // GLE surf complains about `#' symbol
    fprintf(file,   "# central difference update (uncond. unstable)\n"                          ); 
    fprintf(file,   "# N = %d, dt = %.4f, dx = %.4f \n",                              N, dt, dx ); 
    fprintf(file,   "# \n" ); 
    fprintf(file,   "# iteration, coordinate, scalar \n" ); 
*/
    // OUTPUT                                           iteration,    coordinate,       scalar
    for (int i=0; i<N; i++) {    
      fprintf(file,   "%.8f, %.8f, %.8f\n",          ((double) n),  ((double) i),        w[i]   ); 
    }
  } 
  fclose(file);                                                                             return;
}


void eval_x (int Nt, double del_t) { dt = del_t;

  init();
  sprintf(fname, "out/data/t=%.3f, N=%d.dat", ( (double) Nt )*dt, Nt);
  file = fopen(fname, "w+");

  for (int n=0; n<Nt;n++) {     // loop through once

    lax(x, w, u, F); 
    fluxL(x, w, u, F);

  }
  for (int i=0; i<N; i++)       // then write to file

      {    fprintf(file,       "%.8f, %.8f\n", x[i], w[i] );    }

  fclose(file);                                                                             return;
}

/*-----------------------------------------------------------------------------------------------*/

/*int main() {*/

  
  /*eval_x(200, .02);*/
  /*eval_x(400, .02);*/
  /*eval_x(600, .02);*/
  /*eval_x(800, .02);*/
  /*eval_x(1000, .02);*/

  /*[>eval_surf(20);<]*/

  /* */
  /*eval_x(40, .2);*/
  /*eval_x(80, .1);*/
  /*eval_x(160, .05);*/
  /**/

  /*printf("Coordinates:  X = (%.3f, %.3f), \t N=%d", a, b, N);*/
  /*printf("\n (dt, dx) = (%.3f, %.3f)", dt, dx);*/
  /*printf("\n    dt/dx = %.7f \n", dt/dx);*/

  /*return 0;*/
/*}*/


