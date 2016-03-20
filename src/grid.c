/*
 * Author: greg jackson
 * Date: Dec 05 2015
 * (1+1)-dimensional advection
 * with: central diff
 *       upwind
 *       shasta
 */
#include "shasta.h"

double *x, *u, *w;

int coord( int i, bnd X ) {           // NB: only treats {-1, 0, ..., SIZE-1, SIZE}
  switch(X) {
    case P: return (i+SIZE)%SIZE;                  // periodic
    case F: return -(i)/SIZE + (SIZE-i-1)/SIZE + i;   // fixed ends
  }
  return 0;                                  // checked: 16.01.06
}

int pixels() {
  int res = 1;
  for (int i=0; i<DIM; i++) { res*=SIZE; }
  return res;
}

double f0(double x) { 
  return 1./(1.+x*x); 
  /*return 1./(exp( x*2. )+1.)+0.000001;*/
}


void init() {
  /*
   *  Initial Conditions:   x    = { x_1, ... , x_SIZE }   -  uniform grid
   *                        u    = { w_1, ... , u_SIZE }   -  velocity field
   *                        w    = { w_1, ... , w_SIZE }   -  scalar quantity
   */
  int N = pixels();
  x = (double *)malloc( DIM*N*sizeof(double) );
  u = (double *)malloc( DIM*N*sizeof(double) );
  w = (double *)malloc( N*sizeof(double) );
  double x0 = 0.;
  for (int i=0; i<N; i++) {
                            x[i] = x0 + dx*((double) i)   ; 
                            u[i] = 1.                   ;
                            /*w[i] = f0( x[i] - 5. )      ;*/
                            w[i] = f0( x[i] - 5. )      ;
  }
  return;
}

