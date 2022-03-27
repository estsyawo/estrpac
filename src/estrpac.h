#include <R.h>
#include <Rinternals.h>
#include <math.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))
#define signum(x) ((x) >=0.0 ? (1):(-(1)))

// declare function prototypes.
void Kern_Esc(double *z, double *Omg, int *n, int *ncz);
void Kern_DL(double *z, double *Omg, int *n, int *ncz);
