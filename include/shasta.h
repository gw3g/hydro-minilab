#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*-----------------------------------------------------------------------------------------------*/
/*                                                                "global" over ALL source files */

#define N 100
//#define a -30.
//#define b 30.

extern double dt, dx, x0;
extern double x[N], u[N], w[N];

/*-----------------------------------------------------------------------------------------------*/

typedef enum bnd {  // boundary conditions
  P,                // periodic
  F,                // fixed
} bnd; 

int coord( int i, bnd X );            // NB: only treats {-1, 0, ..., N-1, N}

//double f0(double);

void init(double *, double *, double *);
void eval_surf (int, double);
void eval_x (int, double);


// central diff
void centre(double *, double);

// upwind
void upwind(double *, double);

//shasta
void lax(double *,double *, double *, bnd);
void fluxL(double *,double *, double *, bnd);
