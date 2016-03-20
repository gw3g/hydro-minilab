/*
 * Author: greg jackson
 * Date: Jan 17 2016
 * (1+1)-dimensional
 * with: shasta
 */

#include "adv.h"

/*-----------------------------------------------------------------------------------------------*/

/* external parameters */

double dt = .05;
double x[N], u[N], w[N];
double dx = (b-a)/( (double) N );

/*-----------------------------------------------------------------------------------------------*/

double f0(double x) { 
  /*return 1./(1.+x*x); */
  return 1./(exp( x*2. )+1.)+0.000001;
}

double *ev( double t0, double ti ) {
  // get e-den & vel from em tensor
  double s = ti*ti, ee;
  double *res = (double *)malloc( 2*sizeof(double) );
  if (s<1e-10) { res[0] = t0; res[1] = 0.; }
  else {
    ee = sqrt( 4.*t0*t0 - 3.*s ) - t0;
    s /= t0 - ee;
    res[0] = ee; res[1] = ti/s;
  }
  return res;
}

double *tt( double ee, double vv ) {
  // get em tensor from e-den & vel 
  double pp = ee/3., g2 = 1./(1.-vv*vv);
  double *res = (double *)malloc( 2*sizeof(double) );
  res[0] = (ee+pp)*g2 - pp ;
  res[1] = (ee+pp)*g2*vv  ;
  return res;
}

int coord( int i, bnd X ) {                  // NB: only treats {-1, 0, ..., N-1, N}
  switch(X) {
    case P: return (i+N)%N;                  // periodic
    case F: return -(i)/N + (N-i-1)/N + i;   // fixed ends
  }
  return 0;                                  // checked: 16.01.06
}

void ic(double x[], double w[], double u[]) {
  /*
   *  Initial Conditions:   x    = { x_1, ... , x_N }   -  uniform grid
   *                        u    = { w_1, ... , u_N }   -  velocity field
   *                        w    = { w_1, ... , w_N }   -  scalar quantity
   */ 
  double *ev_temp = (double *)malloc( 2*sizeof(double) );;
  for (int i=0; i<N; i++) {
                            x[i] = a + dx*((double) i)   ; 
                            ev_temp = tt( f0( fabs(x[i])-10.), 0. );
                            w[i]   =  ev_temp[0];
                            u[i]   =  ev_temp[1];
                            /*u[i] = 0.                   ;*/
                            /*w[i] = f0( x[i] - 5. )      ;*/
                            /*w[i] = f0( x[i] - 10. )      ;*/
  }
  return;
}

void evo(double x[], double t0[], double tr[], double t) {
  double *ev_temp = (double *)malloc( 2*sizeof(double) );;
  double *boa = (double *)malloc( N*sizeof(double) );
  double *s0  = (double *)malloc( N*sizeof(double) );
  double *sr  = (double *)malloc( N*sizeof(double) );
  double *e   = (double *)malloc( N*sizeof(double) );
  double *v   = (double *)malloc( N*sizeof(double) );
  int ip, im;

  // needs fixing! rather make lax and fluxL return pointers!!
  for (int i=0; i<N; i++) {
    ev_temp = ev( t0[i], tr[i] );
    e[i]   =  ev_temp[0];
    v[i]   =  ev_temp[1];
    boa[i] = tr[i] / t0[i];
  }
  for (int i=0; i<N; i++) {
    ip = coord(i+1,F);
    im = coord(i-1,F);

    s0[i]  = -  x[i]*e[i ] / (3.);
    sr[i]  = -t*x[i]*( e[ip]-e[im] )/(3.*2.*dx);
    t0[i] *= t*x[i];
    tr[i] *= t*x[i];
  }

  lax(x, t0, boa, F);
  lax(x, tr, v, F);

  for (int i=0; i<N; i++) {
    t0[i] += dt*s0[i];
    tr[i] += dt*sr[i];
  }

  fluxL(x, t0, boa, F);
  fluxL(x, tr, v, F);

  for (int i=0; i<N; i++) {
    t0[i] /= t*x[i];
    tr[i] /= t*x[i];
  }

  free(boa);free(s0);free(sr);free(e);free(v);free(ev_temp);
  return;
}

/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[40];

void eval_x (int Nt, double del_t) { dt = del_t;
  double *ev_temp = (double *)malloc( 2*sizeof(double) );
  double *e   = (double *)malloc( N*sizeof(double) );
  double *v   = (double *)malloc( N*sizeof(double) );

  ic(x,w,u);
  sprintf(fname, "out/radial t=%.3f, N=%d.dat", ( (double) Nt )*dt, Nt);
  file = fopen(fname, "w+");

  printf("%.6f, %.6f\n", w[32], u[34]);

  for (int n=1; n<Nt;n++) {
    evo(x,w,u,( (double) n )*dt);
  }

  for (int i=0; i<N; i++) {
    ev_temp = ev( w[i], u[i] );
    e[i]   =  ev_temp[0];
    v[i]   =  ev_temp[1];

    fprintf(file,       "%.8f, %.8f, %.8f\n", x[i], w[i], u[i] );
    /*printf(       "%.8f, %.8f, %.8f\n", x[i], e[i], v[i] );*/
  }

  free(ev_temp);free(e);free(v);

  fclose(file);                                                                             return;
}

void eval_e(int Nt, double del_t, int i) { dt = del_t;
  /*
   * output { time, e-den(x[i]) }
   *
   */
  double *ev_temp = (double *)malloc( 2*sizeof(double) );
  double *e   = (double *)malloc( N*sizeof(double) );
  double *v   = (double *)malloc( N*sizeof(double) );

  ic(x,w,u);
  sprintf(fname, "out/radial i=%d, N=%d.dat", i, Nt);
  file = fopen(fname, "w+");

  for (int n=1; n<Nt;n++) {
    ev_temp = ev( w[i], u[i] );
    e[i]   =  ev_temp[0];
    v[i]   =  ev_temp[1];

    evo(x,w,u,( (double) n )*dt);
    fprintf(file,       "%.8f, %.8f\n", ( (double) n)*dt, e[i] );
  }

  free(ev_temp);free(e);free(v);

  fclose(file);                                                                             return;
}

/*-----------------------------------------------------------------------------------------------*/

int main() {

  eval_x(10, .01);
  eval_x(110, .01);
  eval_x(210, .01);
  eval_x(310, .01);
  eval_x(410, .01);
  eval_e(10000, .01, 202);

  /*eval_x(400, .02);*/
  /*eval_x(600, .02);*/
  /*eval_x(800, .02);*/
  /*eval_x(1000, .02);*/

  /*eval_3d(20);*/

  /* 
  eval_x(40, .2);
  eval_x(80, .1);
  eval_x(160, .05);
  */

  printf("Coordinates:  X = (%.3f, %.3f), \t N=%d", a, b, N);
  printf("\n (dt, dx) = (%.3f, %.3f)", dt, dx);
  printf("\n    dt/dx = %.7f \n", dt/dx);

  return 0;
}


