#include "auto_f2c.h"
#include <math.h>
int func (integer ndim, double *u, const integer *icp, 
     double *par, integer ijac, double *f, double *dfdu, double *dfdp)
{
	/* System generated locals */
    integer dfdu_dim1 = ndim, dfdp_dim1 = ndim;
	/* User defined locals */
	double  x;
	double  gamma;
	x=u[0];
	gamma=par[0];
	f[0] = gamma-x*x;
	return 0;
} 
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, double t, double *u, double *par)
{
	par[0]=0;
	u[0]=0;
	return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int pvls(integer ndim, const doublereal *u,doublereal *par)
{return 0;} /* pvls */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd () { return 0; }
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd () { return 0; } 
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt() { return 0; }
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
