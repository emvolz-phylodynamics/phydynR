/*
 * @author Erik M Volz
 * @date April 18 2014
 * 
 * Equations for lineage states and cumulative hazard of coalescent.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static double *parms; 

//macros
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) < (b)) ? (b) : (a))

#define m (int)parms[0]
#define treeT parms[1]
#define hres parms[2]
#define  Atotal parms[3]

#define F(i, k,l) parms[(int)(4 + hres * (l*m +k ) + i)]
#define fend (int)(4 + hres * pow(m,2))
#define G(i, k,l) parms[fend + (int)(hres * (l*m +k ) + i)]
#define gend (int)(4 + 2 * hres * pow(m,2))
#define Y(i, k)  parms[gend + (int)(hres*k+i)]

/* initializers  */
void initfunc(void (* odeparms)(int *, double *))
{
	int Nparms;
	DL_FUNC get_deSolve_gparms;
	SEXP gparms;
	get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
	gparms = get_deSolve_gparms();
	Nparms = LENGTH(gparms);
	parms = REAL(gparms);
}



/* derivatives */
#define Q(k,l) y[l*m + k]
#define A(k) y[(int)pow(m,2) + k]
#define L y[(int)pow(m,2) + m]
#define dQ(k,l) ydot[l*m + k]
#define dA(k) ydot[(int)pow(m,2) + k]
#define dL ydot[(int)pow(m,2) + m]

void dQAL( int *neq, double *t, double *y, double *ydot, double *yout, int*ip)
{
	// time index
	//~ int i =  (int)min( 1+(int)( hres * (*t) / treeT ), hres);
	int i =  (int)min( (int)( hres * (*t) / treeT ), hres-1);
	int k,l,z,w;
	
	double a[m]; //normalized nlft 
	double sumA = 0.; 
	for (k = 0; k < m; k++) sumA += A(k);
	double r = Atotal / sumA;
	for (k = 0; k < m; k++) { 
		dA(k) = 0.;
		if (Y(i,k) > 0) {
			a[k] = r *  A(k)/ Y(i,k);
			//~ a[k] =  max(0, min(1, r *  A(k)/Y(i,k)));
			//~ a[k] = max( min(r * A(k)/Y(i,k), 1), 0) ;
		} else{
			a[k] = 1.; //
		} 
	}
	
	//dA
	for (k = 0; k < m; k++){
		for (l = 0; l < m; l++){
			if (k==l){
				dA(k) -= a[l] * (F(i,l,k)) * a[k];
			} else{
				dA(k) += ( max(0, (1 - a[k])) * F(i,k,l) + G(i,k,l)) * a[l] ;
				dA(k) -= (F(i,l,k) + G(i,l,k)) * a[k];
			}
		}
	}
	//dQ
	for (z = 0; z < m; z++){ // col of Q
		for (k = 0; k < m; k++){ //row of Q
			dQ(k,z) = 0.;
			for (l = 0. ; l < m; l++){
				if (k!=l){
					if ( Q(l,z) > 0)
					{
						dQ(k,z) += (F(i,k,l) + G(i,k,l)) *  Q(l,z)/  max(Q(l,z), Y(i,l));
					}
					if (Q(k,z) > 0)
					{
						dQ(k,z) -= (F(i,l,k) + G(i,l,k)) *  Q(k,z)/  max(Q(k,z), Y(i,k));
					}
				}
				// coalescent:
				//~ dQ(k,z) -= (F(i,k,l)+F(i,l,k)) * a[l] * Q(k,z)/Y(i,k);
				if (Q(k,z) > 0){
					dQ(k,z) -= F(i,k,l) * a[l] * Q(k,z)/  max(Q(k,z), Y(i,k));
				}
			}
		}
	}
	//dL
	dL = 0.;
	double Ydenom;
	for (k = 0; k < m; k++){
		for (l =0 ; l < m; l++){
			if (k==l){
				Ydenom = max( (Y(i,k)-1),(r*A(k)-1) );
				if (Ydenom > 0)
				{
					dL += max(0,min(1,a[k])) * (  (r*A(k)-1)/Ydenom) * F(i,k,l);
				}
			} else{
				dL += max(0,min(1,a[k])) * max(0,min(1,a[l])) * F(i,k,l);
			}
		}
	}
	dL = max(dL, 0);
}

void dQAL2012( int *neq, double *t, double *y, double *ydot, double *yout, int*ip)
{
	//p_i from 2012 paper
	// time index
	//~ int i =  (int)min( 1+(int)( hres * (*t) / treeT ), hres);
	int i =  (int)min( (int)( hres * (*t) / treeT ), hres-1);
	int k,l,z,w;
	
	double a[m]; //normalized nlft 
	double sumA = 0.; 
	for (k = 0; k < m; k++) sumA += A(k);
	double r = Atotal / sumA;
	for (k = 0; k < m; k++) { 
		dA(k) = 0.;
		if (Y(i,k) > 0) {
			a[k] = r * A(k)/Y(i,k);
			//~ a[k] = max( min(r * A(k)/Y(i,k), 1), 0) ;
		} else{
			a[k] = 1.; //
		} 
	}
	
	//dA
	for (k = 0; k < m; k++){
		for (l = 0; l < m; l++){
			if (k==l){
				dA(k) -= a[l] * (F(i,l,k)) * a[k];
			} else{
				dA(k) += ((1 - a[k]) * F(i,k,l) + G(i,k,l)) * a[l] ;
				dA(k) -= (F(i,l,k) + G(i,l,k)) * a[k];
			}
		}
	}
	//dQ
	for (z = 0; z < m; z++){ // col of Q
		for (k = 0; k < m; k++){ //row of Q
			dQ(k,z) = 0.;
			for (l = 0. ; l < m; l++){
				if (k!=l){
					dQ(k,z) += (F(i,k,l)*(1-a[k]) + G(i,k,l)) * Q(l,z)/Y(i,l);
					dQ(k,z) -= (F(i,l,k)*(1-a[l]) + G(i,l,k)) * Q(k,z)/Y(i,k);
				}
				// coalescent:
				//~ dQ(k,z) -= (F(i,k,l)+F(i,l,k)) * a[l] * Q(k,z)/Y(i,k);
				//dQ(k,z) -= F(i,k,l) * a[l] * Q(k,z)/Y(i,k);
			}
		}
	}
	//dL
	dL = 0.;
	for (k = 0; k < m; k++){
		for (l =0 ; l < m; l++){
			if (k==l){
				dL += a[k] * (  (r*A(k)-1)/(Y(i,k)-1)) * F(i,k,l);
			} else{
				dL += a[k] * a[l] * F(i,k,l);
			}
		}
	}
	
}
