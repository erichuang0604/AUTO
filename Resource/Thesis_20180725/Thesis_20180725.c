#include "auto_f2c.h"
#include <math.h>

int func (integer ndim, double *u, const integer *icp, 
     double *par, integer ijac, double *f, double *dfdu, double *dfdp)
{
	/* System generated locals */
    integer dfdu_dim1 = ndim, dfdp_dim1 = ndim;
	/* User defined locals */
	double  x[4]={0};
	double	Q_1 = 0;
	/* System parameter*/
	double	k_pw = 0.4, k_pv = 0.3, k_qw = -0.03,
			k_qv = -2.8, k_qv2 = 2.1, T = 8.5,
			P_0 = 0.6, Q_0 = 1.3, M = 0.3,
			Y_0 = 20.0, theta_0 = -5.0, E_0 = 1.0,
			C = 12.0, Y_m = 5.0, theta_m = -5.0,
			E_m = 1.0, P_m = 1.0, d_m = 0.05;
			
	double	E_0_prime = 0, Y_0_prime = 0, roll = 0;
	double	P = 0, Q = 0;
	/* Varying parm */
	double	n = 1;
	double	P_1 = 0;
			
	x[0] = u[0];	x[1] = u[1];	x[2] = u[2];	x[3] = u[3];
	Q_1=par[0];
	//Q_1=par[10];
	
	E_0_prime = E_0*(pow((1+pow(C,2)*pow((1/Y_0),2) - 2*C*(1/Y_0)*cos((theta_0/180)*M_PI)),-0.5));
	Y_0_prime = Y_0*(pow((1+pow(C,2)*pow((1/Y_0),2) - 2*C*(1/Y_0)*cos((theta_0/180)*M_PI)),0.5));
	roll = atan((C*(1/Y_0)*sin((theta_0/180)*M_PI))/(1 - C*(1/Y_0)*cos((theta_0/180)*M_PI)));
	
	P = -E_0_prime*Y_0_prime*x[3]*sin(x[2]+(theta_0/180)*M_PI) - (1/n)*E_m*Y_m*x[3]*sin(x[2] - x[0] + (theta_m/180)*M_PI) +
		(Y_0_prime*sin(roll+(theta_0/180)*M_PI) + (1/(pow(n,2)))*Y_m*sin((theta_m/180)*M_PI))*(pow(x[3],2));
	Q = E_0_prime*Y_0_prime*x[3]*cos(x[2]+(theta_0/180)*M_PI) + (1/n)*E_m*Y_m*x[3]*cos(x[2] - x[0] + (theta_m/180)*M_PI) -
		(Y_0_prime*cos(roll+(theta_0/180)*M_PI) + (1/(pow(n,2)))*Y_m*cos((theta_m/180)*M_PI))*(pow(x[3],2));
	
	
	f[0] = x[1];
	f[1] = (P_m - d_m*x[1] + (1/n)*E_m*Y_m*x[3]*sin(x[2] - x[0]- (theta_m/180)*M_PI) + (pow(E_m,2))*Y_m*sin(theta_m/180)*M_PI)/M;
	f[2] = (-k_qv2*(pow(x[3],2)) - k_qv*x[3] + Q - Q_0 - Q_1)/k_qw;
	f[3] = (k_pw*k_qv2*(pow(x[3],2)) + (k_pw*k_qv - k_qw*k_pv)*x[3] + k_qw*(P - P_0 - P_1) - k_pw*(Q - Q_0 - Q_1))/(T*k_qw*k_pv);
	
	return 0;
} 
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, double t, double *u, double *par)
{
	par[0]=1.1322752869E+01;
	u[0]=2.2305977101E-01;	u[1] = 0.0000000000E+00; u[2] = 1.2756854415E-02; u[3] = 9.1866852025E-01;
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