#include "../include/shasta.h"


/*-----------------------------------------------------------------------------------------------*/

/* external parameters */

double dt = .2;
double dx = .4;
/*(b-a)/( (double) N );*/

/*-----------------------------------------------------------------------------------------------*/

int main() {

  /*eval_x(200, .02);*/
  /*eval_x(400, .02);*/
  /*eval_x(600, .02);*/
  /*eval_x(800, .02);*/
  /*eval_x(1000, .02);*/

  /*dt = .2;*/
  eval_surf(20, .4);

  /* 
  eval_x(40, .2);
  eval_x(80, .1);
  eval_x(160, .05);
  */

  printf("Coordinates:  X = (%.3f, %.3f), \t N=%d", x0, x0+N*dx, N);
  printf("\n (dt, dx) = (%.3f, %.3f)", dt, dx);
  printf("\n    dt/dx = %.7f \n", dt/dx);

  return 0;
}

