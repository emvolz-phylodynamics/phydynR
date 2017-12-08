// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>

#include <RcppArmadillo.h>



//~ using namespace arma;
using namespace Rcpp; 
using namespace std; 



const double MIN_Y = 1. ;

typedef std::vector<double> state_type; 

class DPikL0{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
	arma::vec A0; 
	double sumA; 
	int nextant; 
public:
	DPikL0( const List Fs, const List  Gs, const List Ys, const int m, const double hres, const double treeT, const int nextant) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT),nextant(nextant) { };
	
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
    {
		int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); 
		arma::mat F = as<arma::mat>(Fs[i]); 
		arma::mat G = as<arma::mat>(Gs[i]); 
		arma::vec Y = clamp(as<arma::vec>(Ys[i]), MIN_Y, INFINITY ); 
		
		int k,l,z,w;
		
		arma::mat R = (F + G);
		R.each_row() /= Y.t()  ;
		R.diag().zeros(); 
		R.diag() = -arma::sum( R, 0).t(); 
		
		arma::mat Pik = x2Pik( x ); 
		arma::mat dPik = arma::zeros<arma::mat>( m, nextant); 
		double dL = 0.; 
		
		arma::vec A = arma::sum(Pik, 1 );
		arma::vec Ampik = arma::zeros(m); 
		arma::vec Ampik_Y = arma::zeros(m); 
		for (z = 0; z < nextant; z++){
			dPik.col(z) = R * Pik.col(z) ; 
			
			Ampik = arma::clamp( A - Pik.col(z), 0., INFINITY);

			Ampik_Y = arma::clamp( A / Y , 0., 1e6 );
			for (k = 0; k < m; k++){
				for( l = 0; l < m; l++){
					dPik(k,z) -= (Pik(k,z)/Y(k)) * F(k,l) * Ampik_Y(l) ; 
					dL += (Pik(k,z)/Y(k)) * F(k,l) * Ampik_Y(l) ; 
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
	
	arma::mat x2Pik( state_type x ){
		return arma::mat( &x[0], m, nextant ); 
	}
	
	state_type generate_initial_conditions( arma::mat pik0, int m, int nextant ){
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
List solvePikL0(arma::vec times, List Fs, List Gs, List Ys
 , double h0
 , double h1
 , arma::mat pik0
 , double step_size_res
){
	double treeT = std::abs( times(0) - times(times.size()-1)); 
	int m = pik0.n_rows;
	int nextant = pik0.n_cols; 
	double hres = times.size();
	
	DPikL0 dpikl0(Fs, Gs, Ys, m, hres, treeT, nextant ); 
	state_type x = dpikl0.generate_initial_conditions( pik0 , m, nextant); 
	
	using namespace boost::numeric::odeint; 
	if (step_size_res < 1){
		size_t steps = integrate( dpikl0 ,  x , h0 , h1 
				  , (h1-h0)/100. );  
	} else{
		runge_kutta_cash_karp54<state_type> stepper;
		double dt = std::abs( times(1) - times(0)); 
		integrate_const( stepper, dpikl0, x, h0, h1,  std::min( dt/step_size_res, (h1-h0)/step_size_res) );
	}
	arma::mat pik1 = dpikl0.x2Pik( x); 
	double L = x[nextant*m + 1 - 1]; 
	
	List o; 
	o["pik"] = pik1;
	o["L"] = L;
	return o; 
}
