#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*-----------------------------------------------------------------------------------------------*/
/*                                                                "global" over ALL source files */

#define SIZE 100
#define DIM  1

//#define a -30.
//#define b 30.

extern double dt, dx;
extern double *x, *u, *w;;

/*-----------------------------------------------------------------------------------------------*/

typedef enum bnd {  // boundary conditions
  P,                // periodic
  F,                // fixed
} bnd; 

int coord( int i, bnd X );            // NB: only treats {-1, 0, ..., N-1, N}

//double f0(double);

void init();
void eval( int, double);

void centre_update( double *, double *, bnd); // unstable

void upwind_update( double *, double *, bnd); // diffusive

//shasta
void LW_update(     double *, double *, bnd); // diffusive
void flux_correct(  double *, double *, bnd);
void shasta_1d(     double *, double *, bnd);

void evo_tr(double *, double *, double *, double);

// hydro
double *ev( double, double);
double *tt( double, double);
