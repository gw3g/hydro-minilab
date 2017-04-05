/*
 * Author: greg jackson
 * Date: Mar 20 2015
 * radial hydro (t,r) @ midrapidity
 * with: shasta
 */
#include "shasta.h"

/*-----------------------------------------------------------------------------------------------*/
/*                                                          { T00, T0r }  <-->  { ene, vel }     */
double *ev( double t0, double ti ) {
  // get e-den & vel from em tensor
  double s = ti*ti, ee;
  double *res = (double *)malloc( 2*sizeof(double) );
  if (fabs(t0-fabs(ti)) < 1e-6) 
       { res[0] =  3.*(t0-fabs(ti));    res[1] =( 3.*fabs(ti)-2.*t0 )/ti; }
  if (s<1e-6) 
       { res[0] = t0;                   res[1] = .0;                      }
  else { ee = sqrt( 4.*t0*t0 - 3.*s ) - t0;
         res[0] = ee;                   res[1] = (t0-ee)*ti/s;            }
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
  double *boa = (double *)malloc( SIZE*sizeof(double) );
  double *s0  = (double *)malloc( SIZE*sizeof(double) );
  double *sr  = (double *)malloc( SIZE*sizeof(double) );
  double *e   = (double *)malloc( SIZE*sizeof(double) );
  double *v   = (double *)malloc( SIZE*sizeof(double) );
  int ip, im;

  /*printf(" %g, %g, %g \n", x[10], t0[10], tr[10]);*/

  // needs fixing! rather make lax and fluxL return pointers!!
  for (int i=0; i<SIZE; i++) {
    ev_temp = ev( t0[i], tr[i] );
    e[i]   =  ev_temp[0];
    v[i]   =  ev_temp[1];
    boa[i] = tr[i] / t0[i];
  }
  for (int i=0; i<SIZE; i++) {
    ip = coord(i+1,F);
    im = coord(i-1,F);

    s0[i]  = 0.;
    sr[i]  = -t*( e[ip]-e[im] )/(3.*2.*dx);
    /*s0[i]  = -  e[i]/(3.) - t*tr[i]/x[i];*/
    /*sr[i]  = -t*( e[ip]-e[im] )/(3.*2.*dx) - v[i]*t*tr[i]/x[i];*/
    /*t0[i] *= t;*/
    /*tr[i] *= t;*/
  }

  printf(" %g, %g, %g \n", x[50], t0[50], tr[50]);
  LW_update(boa, t0, F);
  LW_update(  v, tr, F);
  printf(" %g, %g, %g \n", x[50], t0[50], tr[50]);

  for (int i=0; i<SIZE; i++) {
    t0[i] += dt*s0[i];
    tr[i] += dt*sr[i];
  }

  flux_correct(boa, t0, F);
  flux_correct(  v, tr, F);

  for (int i=0; i<SIZE; i++) {
    t0[i] /= t;
    tr[i] /= t;
  }

  free(boa);free(s0);free(sr);free(e);free(v);free(ev_temp);
  return;
}

