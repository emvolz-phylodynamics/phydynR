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

class DPikL1{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
	arma::vec A0; 
	double sumA; 
	int nextant; 
public:
	DPikL1( const List Fs, const List  Gs, const List Ys, const int m, const double hres, const double treeT, const int nextant) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT),nextant(nextant) { };
	
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
    {
		int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); 
		arma::mat F = as<arma::mat>(Fs[i]); 
		arma::mat G = as<arma::mat>(Gs[i]); 
		arma::vec Y = clamp(as<arma::vec>(Ys[i]), MIN_Y, INFINITY ); 
		
		int k,l,z,w;
		
		arma::mat R = arma::zeros<arma::mat>(m,m); 
		arma::mat Pik = x2Pik( x ); 
		arma::mat dPik = arma::zeros<arma::mat>( m, nextant); 
		double dL = 0.; 
		
		// normalised P_ik:
		arma::vec colsumPik = arma::sum( Pik, 0).t(); 
		arma::mat nPik = arma::zeros<arma::mat>( m, nextant); 
		for (z = 0; z < Pik.n_cols; z++){
			double s =  sum(Pik.col(z));
			if (s > 0) 
			  nPik.col(z) = Pik.col(z) / s ; 
		}
		
		arma::vec nA = arma::sum(nPik, 1 );
		arma::vec na = nA/Y; 
		
		arma::vec u = arma::zeros(m);
		
		arma::mat phi = arma::repmat(F,1,1); //
		phi.each_row() /= Y.t(); 
		phi.each_col() /= Y; 
		
		for (z = 0; z < nextant; z++){
			//~ u = arma::clamp( 1. - (na - nPik.col(z)/Y) , 0., INFINITY); 
			u = arma::clamp( 1. - na , 0., INFINITY); 
			R = arma::repmat(F,1,1); //
			R.each_col() %= u; 
			R += G; 
			R.each_row() /= Y.t()  ;
			R.diag().zeros(); 
			R.diag() = -arma::sum( R, 0).t(); 
			
			dPik.col(z) = R * Pik.col(z) ; 
			
			for (k = 0; k < m; k++){
				for( l = 0; l < m; l++){
					dPik(k,z) -= Pik(k,z) * (phi(k,l)+phi(l,k)) * (nA(l)- nPik(l,z)); 
					dL += nPik(k,z) * (phi(k,l)+phi(l,k)) * (nA(l)- nPik(l,z))/ 2.;  // note normed pik(k,z) and /2
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
List solvePikL1(arma::vec times, List Fs, List Gs, List Ys
 , double h0
 , double h1
 , arma::mat pik0
 , double step_size_res
){
	double treeT = std::abs( times(0) - times(times.size()-1)); 
	int m = pik0.n_rows;
	int nextant = pik0.n_cols; 
	double hres = times.size();
	
	DPikL1 dpikl0(Fs, Gs, Ys, m, hres, treeT, nextant ); 
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
