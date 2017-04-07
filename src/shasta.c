#include "shasta.h"

void centre_update(double *u, double *w, bnd X) {

  int ic, ip, im;
  for (int i=0; i<SIZE; i++) {
    ic=coord(i  ,X);
    ip=coord(i+1,X);
    im=coord(i-1,X);
    w[ic] -= u[ic]*dt*( w[ip] - w[im] )/(2.*dx);
  }

  return;
}

void upwind_update(double *u, double *w, bnd X) {

  int ic, im;
  for (int i=0; i<SIZE; i++) {
    ic=coord(i  ,X);
    im=coord(i-1,X);
    w[ic] -= u[ic]*dt*( w[ic] - w[im] )/(dx);
  }

  return;
}


void LW_update(double *u, double *w, bnd X) {
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

  for (int i=-1; i<SIZE; i++) {
    Wc = Wn; Wn=0.;                                       // assign Wc > Wn and reset Wn

    ic=coord(i  ,X);
    ip=coord(i+1,X);

    Qp = ( .5 - l*u[ic] )/( 1. + l*( u[ip] - u[ic] ) );
    Qm = ( .5 + l*u[ip] )/( 1. + l*( u[ip] - u[ic] ) );    // evaluated at i+1

                                    Di = w[ip] - w[ic];    // BCs here

    if (i>-1)    {Wc += Qp*( .5*Qp*Di + w[ic] ); w[i]=Wc;} // if ... add to   Wc 
    if (i<+SIZE) {Wn -= Qm*( .5*Qm*Di - w[ip] );         } // if ... sub from Wn

  }
                                                           return;
}

void flux_correct(double *u, double *ws, bnd X) {  
  /*
   * flux correct :     ws = diffusive soln.
   *                    u  = velocity field
   *                    X  = boundary cond.
   *
   */
  double                  Ai;                               // anti-diffusive flux
  int        im,    ic,   ip;                               // index triple <i-1,i,i+1>
  int                    sig;                               // sig_i = +/- 1

  double D[SIZE];                                              // no easy way to avoid? :-(

  for (int i=0; i<SIZE; i++) {
    ip=coord(i+1,X);

    D[i] = ws[ip] - ws[ i ];
  }

  for (int i=0; i<SIZE-1; i++) {
    ic=coord(i  ,X);
    ip=coord(i+1,X); 
    im=coord(i-1,X);

    sig = (int) ( (D[ic]>0) - (D[ic]<0) );                  //  = sign( D[ic] )
    Ai = sig*fmax(0,
                      fmin( fabs(D[ic]*.11), 
                            fmin( sig*D[ip], sig*D[im] ) 
                          )
                    );
    ws[ic]-= Ai;
    ws[ip]+= Ai;
  }
                                                            return;
}

void shasta_1d( double *u, double *w, bnd X) {
  LW_update( w, u, X);
  flux_correct( w, u, X);
}
