/*
 * @author Erik M Volz
 * @date October 14 2014
 * 
 * Equations for lineage through time and distribution of node heights
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

#define pend 3
#define F(i, k,l) parms[(int)(pend + hres * (l*m +k ) + i)]
#define fend (int)(pend + hres * pow(m,2))
#define G(i, k,l) parms[fend + (int)(hres * (l*m +k ) + i)]
#define gend (int)(pend + 2 * hres * pow(m,2))
#define Y(i, k)  parms[gend + (int)(hres*k+i)]
#define yend (int)(gend+hres*m)
#define notSampledYet(i,k) parms[yend + (int)(hres*k+i) ]

/* initializers  */
void initfunc_dCA(void (* odeparms)(int *, double *))
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
#define A(k) y[k]
#define dA(k) ydot[k]
#define dLambda ydot[m]
#define Lambda y[m] 

//~ R version: 
	//~ dA <- function(h, A, parms, ...)
	//~ {
		//~ nsy <- not.sampled.yet(h) 
		//~ with(get.fgy(h), 
		//~ { 
			//~ A_Y 	<- (A-nsy) / .Y
			//~ csFpG 	<- colSums( .F + .G )
			//~ list( .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) )
		//~ })
	//~ }
void dCA( int *neq, double *t, double *y, double *ydot, double *yout, int*ip)
{
	int i =  (int)min( (int)( hres * (*t) / treeT ), hres-1);
	//~ if ((*t) < 10)
	//~ {
		//~ printf( "Fs %f %f \n", F(i, 3, 0), F(i, 0, 3)); 
		//~ printf( "Ys %f %f \n", Y(i, 0), Y(i, 1)); 
		//~ printf( "NSY %f %f \n", notSampledYet(i, 0), notSampledYet(i, 3)); 
		//~ printf( "t %f  treeT %f hres %f \n", (*t), treeT, hres); 
	//~ }
	int k,l;
	
	
	double csFpG[m]; 
	for (k = 0; k < m; k++)
	{
		csFpG[k] = 0.;
	}
	double a[m]; //normalized nlft 
	for (k = 0; k < m; k++) { 
		dA(k) = 0.;
		if (Y(i,k) > 0) {
			a[k] =  max(0, (A(k)-notSampledYet(i,k))) / Y(i,k); //max(1, Y(i,k)) );
		} else{
			a[k] = 1.; //
		} 
		for ( l = 0; l < m; l++)
		{
			csFpG[l] += G(i, k, l) + F(i, k, l);
		}
	}
	
	//dA
	//~ list( .G %*% A_Y - csFpG * A_Y + (.F %*% A_Y) * pmax(1-A_Y, 0) )
	for (k = 0; k < m; k++){
		dA(k) -= a[k] * csFpG[k];
		for (l = 0; l < m; l++){
			dA(k) += a[l] * G(i,k,l);
			dA(k) += max(0,(1-a[k])) * F(i,k,l) * a[l];
			dLambda += F(i,k,l) * a[k] * a[l]; //note- would be a little more accurate to have (A-1)/(Y-1) if k==l
		}
	}
	
	
	
	/*for (k = 0; k < m; k++){
		for (l = 0; l < m; l++){
			if (k==l){
				dA(k) -= a[l] * (F(i,l,k)) * a[k];
			} else{
				dA(k) += ( F(i,k,l) + G(i,k,l)) * a[l] ; //(1 - a[k]) * ??
				dA(k) -= ( F(i,l,k) + G(i,l,k)) * a[k]; //TODO 
			}
		}
	}
	*/
	
	//~ if ((*t) < 10)
	//~ {
		//~ printf("TIME %f\n", (*t)); 
		//~ for (k = 0; k < m; k++)
		//~ {
			//~ printf("%f %f %f %f %f \n", Y(i, k), notSampledYet(i,k), A(k) , a[k], dA(k));
		//~ }
		//~ printf("~~~~~~~~~~~\n");
		//~ for (k = 0; k < m; k++)
		//~ {
			//~ for (l = 0; l < m; l++)
			//~ {
				//~ printf("%f ", F(i,k,l) );
			//~ }
			//~ printf("\n"); 
		//~ }
		//~ printf("~~~~~~~~~~~\n");
		//~ for (k = 0; k < m; k++)
		//~ {
			//~ for (l = 0; l < m; l++)
			//~ {
				//~ printf("%f ", G(i,k,l) );
			//~ }
			//~ printf("\n"); 
		//~ }
		//~ printf("~~~~~~~~~~~\n");
		//~ printf("~~~~~~~~~~~\n");
	//~ }
}
