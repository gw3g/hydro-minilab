#include "shasta.h"


/*-----------------------------------------------------------------------------------------------*/

/* external parameters */

double dt = .2;
double dx = .4;
/*(b-a)/( (double) SIZE );*/

/*-----------------------------------------------------------------------------------------------*/
void eval (int, double);
void ener (int, double);
void init (double *, double *, double *);
/*-----------------------------------------------------------------------------------------------*/

int main() {

  /*eval_x(200, .02);*/
  /*eval_x(400, .02);*/
  /*eval_x(600, .02);*/
  /*eval_x(800, .02);*/
  /*eval_x(1000, .02);*/

  /*dt = .2;*/
  eval(100, .1);
  ener(100, .1);

  /* 
  eval_x(40, .2);
  eval_x(80, .1);
  eval_x(160, .05);
  */

  printf("Coordinates:  X = (%.3f, %.3f), \t SIZE=%d", x[0], x[SIZE-1], SIZE);
  printf("\n (dt, dx) = (%.3f, %.3f)", dt, dx);
  printf("\n    dt/dx = %.7f \n", dt/dx);

  return 0;
}

/*-----------------------------------------------------------------------------------------------*/

FILE *file; char fname[40];

void eval (int Nt, double del_t) { dt = del_t;

  init();
  /*sprintf(fname, "out/data/shasta, tf=%g.dat", ( (double) Nt )*dt);*/
  sprintf(fname, "out/data/Evo.dat", ( (double) Nt )*dt);
  file = fopen(fname, "w+");
  // GLE surf complains about `#' symbol
  fprintf(file,   "# SIZE = %d, dt = %.4f, dx = %.4f \n",                         SIZE, dt, dx ); 
  fprintf(file,   "# \n" ); 
  fprintf(file,   "# t, x, w \n" ); 

  for (int n=0; n<Nt;n++) { 

    /*shasta_1d( u, w, P);       // central diff update*/
    /*LW_update( u, w, P);*/
    /*flux_correct( u, w, P);*/
    evo_tr(x,w,u,.6+dt*((double) n));

    // OUTPUT                                      iteration,    coordinate,       scalar
    for (int i=0; i<SIZE; i++) {    
      fprintf(file,   "%.8f, %.8f, %.8f, %8f\n",       dt*((double) n),  x[i],        w[i], u[i]   ); 
    }
  } 
  fclose(file);                                                                             return;
}

void ener (int Nt, double del_t) { dt = del_t;

  int ix = 50;
  init();
  sprintf(fname, "out/data/e(%g), tf=%g.dat", x[ix], ( (double) Nt )*dt);
  file = fopen(fname, "w+");
  // GLE surf complains about `#' symbol
  fprintf(file,   "# SIZE = %d, dt = %.4f, dx = %.4f \n",                         SIZE, dt, dx ); 
  fprintf(file,   "# \n" ); 
  fprintf(file,   "# t, x, w \n" ); 
  double *temp_ev;

  for (int n=0; n<Nt;n++) { 

    /*shasta_1d( u, w, P);       // central diff update*/
    /*LW_update( u, w, P);*/
    /*flux_correct( u, w, P);*/
    evo_tr(x,w,u,.6+dt*((double) n));

    // OUTPUT                                      iteration,    coordinate,       scalar
    for (int i=0; i<SIZE; i++) {    
      temp_ev = ev( w[ix], u[ix] );
      fprintf(file,   "%.8f, %.8f \n",   .6+dt*((double) n), temp_ev[0]   ); 
    }
  } 
  fclose(file);                                                                             return;
}


