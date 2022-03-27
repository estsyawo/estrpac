#include "estrpac.h"
#include <math.h>

//****************************************************************************//
// compute the n x n kernel for the ICM estimator of Escanciano 2006
void Kern_Esc(double *z, double *Omg, int *n, int *ncz)
{
  double cp,Dij,Dijsq,Djl,Djlsq,Dil,Dilsq,dd,Aijl;
  int i,j,l,k;
    cp = (double) pow(M_PI,((*ncz/2) - 1))/tgamma((*ncz/2) + 1) ;
  for( i=0; i<*n; i++){
    for ( j=0; j<=i; j++) {
        Aijl=0.0;
      for ( l=0; l<*n; l++) {
          Dijsq=0.0;Djlsq=0.0;Dilsq=0.0; dd=0.0;
          for (k=0; k<*ncz; k++) {
              Dij=0.0;Djl=0.0;Dil=0.0;
              Dij = z[k*(*n)+i]-z[k*(*n)+j]; Dijsq += Dij*Dij;
              Djl = z[k*(*n)+j]-z[k*(*n)+l]; Djlsq += Djl*Djl;
              Dil = z[k*(*n)+i]-z[k*(*n)+l]; Dilsq += Dil*Dil;
              dd += Dil*Djl;
          }//end for (k=0; k<*ncz; k++)
          if ( (Dijsq==0 && Dilsq>0.0) || (Dijsq>0.0 && Dilsq==0.0) || (Dijsq>0.0 && Djlsq==0.0)) {
              Aijl += M_PI;
          }else if(Dijsq==0.0 && Dilsq==0.0){
              Aijl += 2.0*M_PI;
          }else{
              dd = dd/sqrt(Dilsq*Djlsq);
              if (absval(dd)>1.0) { //avoid numerical overflow beyond the [-1,1] bound
                  dd = (double) signum(dd);
              }
              Aijl += absval(M_PI - acos(dd));
          }
      }// end for ( l=0; l<*n; l++)
      Omg[ j*(*n)+i] = cp*Aijl/(*n);
      if (j!=i) {
          Omg[i*(*n)+j]=Omg[j*(*n)+i];//by symmetry
      }//end if(j!=i)
    }// end for ( j=0; j< i; j++)
  }// for( i=1; i< *n;i++)
}// end function
//****************************************************************************//


//****************************************************************************//
// compute the n x n kernel for the ICM estimator of Dominguez & Lobato 2004
void Kern_DL(double *z, double *Omg, int *n, int *ncz)
{
  double tmp,Aijl;
  int i,j,l,k;
  for( i=0; i<*n; i++){
    for ( j=0; j<=i; j++) {
      Aijl=0.0;
      for ( l=0; l<*n; l++) {
        tmp=1.0;// initialise to 1
        for (k=0; k<*ncz; k++) {
          if((z[k*(*n)+l] < z[k*(*n)+i]) || (z[k*(*n)+l] < z[k*(*n)+j])){
            tmp=0.0;
            break;
          }// end if
        }//end for (k=0; k<*ncz; k++)
          Aijl += tmp;
      }// end for ( l=0; l<*n; l++)
      Omg[ j*(*n)+i] = Aijl/(*n);
      if (j!=i) {
        Omg[i*(*n)+j]=Omg[j*(*n)+i];//by symmetry
      }//end if(j!=i)
    }// end for ( j=0; j<= i; j++)
  }// for( i=1; i< *n;i++)
}// end function
//****************************************************************************//