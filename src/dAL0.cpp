/* lineages through time and hazard of co
 */
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]


static const double A_EPS = 0.00;


using namespace std;
using namespace Rcpp; 


#define A_x(k) std::max(A_EPS, y[ k])

//[[Rcpp::export]]
List dAL(double t, NumericVector y,  List parms){
	//y : (A, L )
	// t : h
	// time index.. note times decreasing
	//~ int it = (int)parms["times.size"] - 1 - (int)floor( t / (double)parms["deltah"] );
	NumericVector times  = as<NumericVector>(parms["times"]);  //TODO this is slow; figure out how to get rid of this
	int tfgylen = times.size();
	int it = std::min(tfgylen-1 , std::max(0, (int)floor( t / (double)parms["deltah"] )));
	int m = (int)parms["m"];
	NumericMatrix F = as<NumericMatrix>(as<List>(parms["Fs"])[it]);
	NumericMatrix G = as<NumericMatrix>(as<List>(parms["Gs"])[it]);
	NumericVector Y = as<NumericVector>(as<List>(parms["Ys"])[it]);
	
	int k,l,z,w;
	
	NumericVector dAL(m+1);
	double a[m]; //normalized nlft 
	double sumA = 0.; 
	for (k = 0; k < m; k++) sumA += A_x(k);
	for (k = 0; k < m; k++) { 
		dAL(k) = 0.;
		if (Y(k) > 0) {
			//~ a[k] = std::max(0., std::min(1., A_x(k)/ Y(k)));
			a[k] = std::max(0., ( A_x(k)/ Y(k)));
		} else{
			//~ a[k] = std::min(1., std::max(0., A_x(k))); 
			a[k] = std::max(0., A_x(k)); 
		} 
	}
	//dA
	//  does not conserve A
	if (sumA > 0 ){
		for (k = 0; k < m; k++){
			for (l = 0; l < m; l++){
				if (k==l){
					dAL(k) -= a[l] * (F(l,k)) * a[k]; //max(0., A_x(k)-1.) / max(1., (Y(k)-1.)); //a[k];
				} else
				{
					//~ dAL(k) += ( std::max(0.0, (1.0 - a[k])) * F(k,l) + G(k,l)) * a[l] ;
					dAL(k) += ( (1.0 - a[k]) * F(k,l) + G(k,l)) * a[l] ;
					//~ dAL(k) += (  F(k,l) + G(k,l)) * a[l] ;
					dAL(k) -= ( F(l,k) + G(l,k)) * a[k];
				}
				dAL(m) += a[l] * (F(l,k)) * a[k];
			}
		}
	}
	return List::create( dAL) ;
}
