/*
 * Author: greg jackson
 * Date: Dec 05 2015
 * (1+1)-dimensional advection
 * with: central diff
 *       upwind
 *       shasta
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*-----------------------------------------------------------------------------------------------*/

/* external parameters */

#define N 100
#define a 0.
#define b 40.

double dt = .05;
double x[N], u[N], w[N];
double dx = (b-a)/( (double) N );

typedef enum bnd {  // boundary conditions
  P,                // periodic
  F,                // fixed
} bnd; 

/*-----------------------------------------------------------------------------------------------*/

double f0(double x) { 
  /*return 1./(1.+x*x); */
  return 1./(exp(x*2.)+1.);
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
void lw(double w[], double u[], bnd X) {
  /*
   * Low order temp solution (Lax-Wrendoff) : w = { w^\bullet_i }
   *                                          u = { u_i }
   *                                          X = boundary cond.
   */
  double l = dt/dx;
                                                          // Dynamic programming approach
  int        ic,    ip;                                   // index pair <i,i+1>
  double     Qp,    Qm,   Di;                             // Q_\pm, \Dela_i
  double     Wc,    Wn = 0.;                              // w {current,next}

  for (int i=-1; i<N; i++) {
    Wc = Wn; Wn=0.;                                       // assign Wc > Wn and reset Wn

    ic=coord(i  ,X);
    ip=coord(i+1,X);

    Qp = ( .5 - l*u[ic] )/( 1. + l*( u[ip] - u[ic] ) );
    Qm = ( .5 + l*u[ip] )/( 1. + l*( u[ip] - u[ic] ) );    // evaluated at i+1

                                    Di = w[ip] - w[ic];    // BCs here

    if (i>-1) {Wc += Qp*( .5*Qp*Di + w[ic] ); w[i]=Wc;}    // if ... add to   Wc 
    if (i<+N) {Wn -= Qm*( .5*Qm*Di - w[ip] );         }    // if ..  sub from Wn

  }
                                                           return;
}

void fc(double ws[], double u[], bnd X) {  
  /*
   * flux correct :     ws = diffusive soln.
   *                    u  = velocity field
   *                    X  = boundary cond.
   *
   */
  double                  Ai;                               // anti-diffusive flux
  int        im,    ic,   ip;                               // index triple <i-1,i,i+1>
  int                    sig;                               // sig_i = +/- 1

  double D[N];                                              // no easy way to avoid? :-(

  for (int i=0; i<N; i++) {
    ip=coord(i+1,X);

    D[i] = ws[ip] - ws[ i ];
  }

  for (int i=0; i<N-1; i++) {
    ic=coord(i  ,X);
    ip=coord(i+1,X); 
    im=coord(i-1,X);

    sig = (int) ( (D[ic]>0) - (D[ic]<0) );                  //  = sign( D[ic] )
    Ai = sig*fmax(0,
                      fmin( fabs(D[ic]/8), 
                            fmin( sig*D[ip], sig*D[im] ) 
                          )
                    );
    ws[ic]-= Ai;
    ws[ip]+= Ai;
  }
                                                            return;
}

void ic() {
  /*
   *  Initial Conditions:   x    = { x_1, ... , x_N }   -  uniform grid
   *                        u    = { w_1, ... , u_N }   -  velocity field
   *                        w    = { w_1, ... , w_N }   -  scalar quantity
   */
  for (int i=0; i<N; i++) {
                            x[i] = a + dx*((double) i)   ; 
                            u[i] = 1.                   ;
                            w[i] = f0( x[i] - 10. )      ;
  }
  return;
}

/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[40];

void eval_3d (int Nt) {

  ic();
  sprintf(fname, "out/3D, t=%.3f-(forward).dat", ( (double) Nt )*dt);
  file = fopen(fname, "w+");

  for (int n=0; n<Nt;n++) { 

    centre( w, 1.);       // central diff update

    // OUTPUT                                           iteration,    coordinate,       scalar
    for (int i=0; i<N; i++) {    
      fprintf(file,   "%.8f, %.8f, %.8f\n",       ((double) n)*dt,          x[i],         w[i]   ); 
    }
  } 
  fclose(file);                                                                             return;
}


void eval_x (int Nt, double del_t) { dt = del_t;

  ic();
  sprintf(fname, "out/t=%.3f, N=%d.dat", ( (double) Nt )*dt, Nt);
  file = fopen(fname, "w+");

  for (int n=0; n<Nt;n++) {     // loop through once

    lw( w, u, F); 
    fc( w, u, F);

  }
  for (int i=0; i<N; i++)       // then write to file

      {    fprintf(file,       "%.8f, %.8f\n", x[i], w[i] );    }

  fclose(file);                                                                             return;
}

/*-----------------------------------------------------------------------------------------------*/

int main() {

  eval_x(200, .02);
  eval_x(400, .02);
  eval_x(600, .02);
  eval_x(800, .02);
  eval_x(1000, .02);

  printf("Coordinates:  X = (%.3f, %.3f), \t N=%d", a, b, N);
  printf("\n (dt, dx) = (%.3f, %.3f)", dt, dx);
  printf("\n    dt/dx = %.7f \n", dt/dx);

  return 0;
}


