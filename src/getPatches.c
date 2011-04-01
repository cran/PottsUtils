#include <R.h>
#include "Rinternals.h"

SEXP getPatches(SEXP sa, SEXP sb, SEXP snbond, SEXP snvert)
{
  int *a = INTEGER(sa); /* transfer a vector*/
  int *b = INTEGER(sb);
  int n1 = asInteger(snbond); /* transfer a scalor*/
  int n2 = asInteger(snvert);
  
  SEXP val;
  PROTECT( val = allocVector(INTSXP, n2)); /* create a new vector for return value*/
  int *f = INTEGER(val);  

  int i, p0, q0, p1, q1;
  
  for(i = 0; i < n2; i++){
    f[i] = i;
  }

  for(i = 0; i < n1; i++){
    p0 = a[i];
    q0 = b[i];
    p1 = f[p0];
    q1 = f[q0];
    while(p1 != q1){
      if(q1 < p1){
	f[p0] = q1;
	p0 = p1;
	p1 = f[p1];
      }
      else{
	f[q0] = p1;
	q0 = q1;
	q1 = f[q1];
      }

    }
  }
  for(i = 0; i < n2; i++){
    f[i] = f[f[i]];
  }

  UNPROTECT(1);
  return val;
}


