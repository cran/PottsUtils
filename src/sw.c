#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

SEXP getPatches(SEXP sa, SEXP sb, SEXP snbond, SEXP snvert)
{
	int *a = INTEGER(sa); 
	int *b = INTEGER(sb);
	int n1 = asInteger(snbond); 
	int n2 = asInteger(snvert);
	

	SEXP val, vertex;
	PROTECT( val = allocVector(INTSXP, n2)); 
	PROTECT( vertex = allocVector(INTSXP, n2)); 
	int *f = INTEGER(val);  
	int *v = INTEGER(vertex);  
	
	int i, p0, q0, p1, q1;	
	
	for(i = 0; i < n2; i++){
		f[i] = i;
		v[i] = i;
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
	
	SEXP patches;
	PROTECT(patches = eval(lang3(install("split"), vertex, val),R_BaseEnv));
		
	UNPROTECT(3);

	return patches;
}


SEXP sw(SEXP sbondProbs, SEXP soneIteration, SEXP sedges, SEXP snedge, 
        SEXP sniter, SEXP snvert, SEXP sncolor)
{
	if (TYPEOF(sbondProbs) != REALSXP)
		error("'bondProbs' must be of type 'double'.");
	if (TYPEOF(soneIteration) != INTSXP)
		error("'oneIteration' must be of type 'integer'.");
	if (TYPEOF(sedges) != INTSXP)
		error("'edges' must be of type 'integer'.");
    if (TYPEOF(snedge) != INTSXP)
		error("'nedge' must be of type 'integer'.");
	if (TYPEOF(sniter) != INTSXP)
		error("'niter' must be of type 'integer'.");
	if (TYPEOF(snvert) != INTSXP)
		error("'nvert' must be of type 'integer'.");
	if (TYPEOF(sncolor) != INTSXP)
		error("'ncolor' must be of type 'integer'.");
     
	
	double *bondProbs = REAL(sbondProbs);
	int *oneIteration = INTEGER(soneIteration);
    int *edges = INTEGER(sedges);
    int nedge = asInteger(snedge);
	int niter = asInteger(sniter);
	int nvert = asInteger(snvert);
	int ncolor = asInteger(sncolor);
	
	int i, j;

	int *bondsPick = (int *) R_alloc(nedge, sizeof(int));

	SEXP scolors = PROTECT(allocMatrix(INTSXP, nvert, niter));
	int *colors = INTEGER(scolors);

	SEXP snbond = PROTECT(allocVector(INTSXP, 1));
	int *nbond = INTEGER(snbond);

	GetRNGstate();
	
	for (i = 0; i < niter; i++) { 

		nbond[0] = 0;

		/* build bonds */
		for (j = 0; j < nedge; j++) {
			double sunif; 
			if (oneIteration[edges[j]] == oneIteration[edges[j+nedge]]){
				sunif = unif_rand();
				if (sunif < bondProbs[j]){
					bondsPick[j] = 1;
					nbond[0] = nbond[0] + 1;
				}
				else{
					bondsPick[j] = 0;
				}
				
			}
			else{
				bondsPick[j] = 0;
			}
		}
		
		if (nbond[0] > 0){
			/* obtain patches */
			SEXP sa = PROTECT(allocVector(INTSXP, nbond[0]));
			SEXP sb = PROTECT(allocVector(INTSXP, nbond[0]));
			int *a = INTEGER(sa);
			int *b = INTEGER(sb);
			int k = 0;
			for (j = 0; j < nedge; j++) {
				if (bondsPick[j] == 1){
					a[k] = edges[j];
					b[k] = edges[j + nedge];
					k = k + 1;
				}
			}
			
			
			SEXP patches = PROTECT(getPatches(sa, sb, snbond, snvert));
			
			/* obtain new colors of each patch */
			//int *newColors = (int *) Calloc(LENGTH(patches), sizeof(int));
			int *newColors = Calloc(LENGTH(patches), int);
			double crand; 
			for( j = 0; j < LENGTH(patches); j++){
				crand = rand() % ncolor;
				newColors[j] = (int) crand; 
			}
			/* assign new colors to each patch */
			for (j = 0; j < LENGTH(patches); j++) {
				int *ver = INTEGER(VECTOR_ELT(patches, j));
				int nver = LENGTH(VECTOR_ELT(patches, j));
				for (k = 0; k < nver; k++) {
					oneIteration[ver[k]] = newColors[j];
					colors[ver[k] + i*nvert] = newColors[j];
				}
			}
		
			Free(newColors);
			UNPROTECT(3);
		}
		else{
			//int *newColors = (int *) Calloc(nvert, sizeof(int));
			int *newColors = Calloc(nvert, int);
			double crand; 
			for( j = 0; j < nvert; j++){
				crand = rand() % ncolor;
				newColors[j] = (int) crand;
			}
		
			for (j = 0; j < nvert; j++) {
				oneIteration[j] = newColors[j];
				colors[j + i*nvert] = newColors[j];
			}
			Free(newColors);
		}
	}

	PutRNGstate();
	
	UNPROTECT(2);

	return scolors;
}

