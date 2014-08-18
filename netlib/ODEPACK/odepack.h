/*
 * This file has been extracted and edited from the ARPACK++ distribution
 * to define the Fortran routines for use in C. The obvious conversions
 * have been used, in agreement with the CBLAS mappings.
 */

// double precision routines.
typedef (*dode_f)(int *neq, double *t, double *y, double *ydot);
typedef (*djac_f)(int *neq, double *t, double *y,
		  int *ml, int *mu, double *pd, int *nrowpd);

void dlsode_(dode_f f, int *neq, double *y, double *t, double *tout,
	     int *itol, double *rtol, double *atol, int *atask,
	     int *istate, int *opt, double *rwork, int *lrw,
	     int *iwork, int *liw, djac_f jac, int *mf);
