// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>

#include <RcppArmadillo.h>



using namespace arma;
using namespace Rcpp; 
using namespace std; 


const double MIN_Y = 1e-12 ;

typedef std::vector<double> state_type; 

class DPikL0{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
	vec A0; 
	double sumA; 
	int nextant; 
public:
	DPikL0( const List Fs, const List  Gs, const List Ys, const int m, const double hres, const double treeT, const int nextant) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT),nextant(nextant) { };
	
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
    {
		int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); 
		mat F = as<mat>(Fs[i]); 
		mat G = as<mat>(Gs[i]); 
		vec Y = clamp(as<vec>(Ys[i]), MIN_Y, INFINITY ); 
		
		int k,l,z,w;
		
		mat R = (F + G);
		R.each_row() /= Y.t()  ;
		R.diag().zeros(); 
		R.diag() = -sum( R, 0).t(); 
		
		mat Pik = x2Pik( x ); 
		mat dPik = zeros<mat>( m, nextant); 
		double dL = 0.; 
		
		vec A = sum(Pik, 1 );
		vec Ampik = zeros(m); 
		for (z = 0; z < nextant; z++){
			dPik.col(z) = R * Pik.col(z) ; 
			
			Ampik = clamp( A - Pik.col(z), 0., INFINITY);
			for (k = 0; k < m; k++){
				for( l = 0; l < m; l++){
					dPik(k,z) -= (Pik(k,z)/Y(k)) * F(k,l) * Ampik(l) / Y(l) ; 
					dL += (Pik(k,z)/Y(k)) * F(k,l) * Ampik(l) / Y(l) ; 
				}
			}
		}
		
		w = 0; 
		for (int z = 0; z < nextant; z++){
			for (int k = 0; k < m; k++){
				dxdt[w] = dPik( k, z); 
				w++;
			}
		}
		dxdt[m*nextant+1-1] = dL; 
	}
	
	mat x2Pik( state_type x ){
		return mat( &x[0], m, nextant ); 
	}
	
	state_type generate_initial_conditions( mat pik0, int m, int nextant ){
		state_type x( 1 + m * nextant , 0.); 
		int w = 0; 
		for (int z =0; z < nextant; z++){
			for (int k = 0; k < m; k++){
				x[w] = pik0(k,z);
				w++; 
			}
		}
		return x;
	}
};

//[[Rcpp::export()]]
List solvePikL0(vec times, List Fs, List Gs, List Ys
 , double h0
 , double h1
 , mat pik0
){
	double treeT = std::abs( times(0) - times(times.size()-1)); 
	int m = pik0.n_rows;
	int nextant = pik0.n_cols; 
	double hres = times.size();
	
	DPikL0 dpikl0(Fs, Gs, Ys, m, hres, treeT, nextant ); 
	state_type x = dpikl0.generate_initial_conditions( pik0 , m, nextant); 
	
	size_t steps = boost::numeric::odeint::integrate( dpikl0 ,  x , h0 , h1 
			  , (h1-h0)/100. );  
	
	
	mat pik1 = dpikl0.x2Pik( x); 
	double L = x[nextant*m + 1 - 1]; 
	
	List o; 
	o["pik"] = pik1;
	o["L"] = L;
	return o; 
}
