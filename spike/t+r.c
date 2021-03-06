/*
 * Author: greg jackson
 * Date: Mar 20 2015
 * radial hydro (t,r) @ midrapidity
 * with: shasta
 */
#include "shasta.h"

double *x, *u, *w;

/*-----------------------------------------------------------------------------------------------*/
/*                                                          { T00, T0r }  <-->  { ene, vel }     */
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

/*-----------------------------------------------------------------------------------------------*/

void evo_tr(double x[], double t0[], double tr[], double t) {
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

