/*
 * Equations on logit-transform of state vectors p_ik
 */

// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>

#include <RcppArmadillo.h>



//~ using namespace arma;
using namespace Rcpp; 
using namespace std; 


const double MIN_Y = 1.;
const double MIN_P = 1e-12;
const double MAX_P = 1.-MIN_P;

typedef std::vector<double> state_type; 



// this version does not solve extra eqns for A; also conserves sum A_k
class DQAL3{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
	arma::vec A0; 
	double sumA; 
public:
	DQAL3( List Fs, List  Gs, List Ys, int m, double hres, double treeT, arma::vec _xA0 ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT),A0(_xA0) {
		sumA = sum(_xA0); 
	};
	
	
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
    {
		// time index
		//~ int i =  (int)min( 1+(int)( hres * (*t) / treeT ), hres);
		//~ NOTE hres = length(times)
		//~ NOTE treeT = max(times)
		int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); 
		arma::mat F = as<arma::mat>(Fs[i]); 
		arma::mat G = as<arma::mat>(Gs[i]); 
		arma::vec Y = arma::clamp(as<arma::vec>(Ys[i]), MIN_Y, INFINITY ); 
		
		int k,l,z,w;
		double pdot; 
		
		arma::vec A = A_from_state( x ) ;
		//~ vec a = clamp( A / Y , 0., 1. ); 
		arma::vec a = A / Y ; 
		
		//d logodds(Q)
		for (z = 0; z < m; z++){ // col of Q
			for (k = 0; k < m; k++){ //row of Q
				pdot = 0.; 
				for (l = 0. ; l < m; l++){
					if (k!=l){
						pdot += (F(k,l) + G(k,l)) *  Q(x, l,z)/  std::max(Q(x,l,z), Y(l)) - (F(l,k) + G(l,k)) *  Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
					}
					// coalescent:
					if (Q(x, k,z) > 0){ 
						pdot -= F(k,l) * a(l) * Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
					}
				}
				dxdt[ Qind(k,z) ] = pdot / ( Q(x,k,z) * (1.-Q(x,k,z) )); 
			}
		}
		//dL
		double dL = 0.;
		double Ydenom;
		for (k = 0; k < m; k++){
			for (l =0 ; l < m; l++){				
				if (k == l && A(k) >= 1. ){
					dL += (A(k) / Y(k)) * ((A(k)-1.) / Y(k)) * F(k,l) ; 
				} else {
					dL += a(k) * a(l) * F(k,l);
				}
			}
		}
//~ cout << A << endl;
//~ cout << A0 << endl;
//~ cout << Y << endl;
//~ cout << a << endl; 
//~ cout << dL << endl; 
//~ cout << x[Lind()] << endl; 
//~ cout << endl; 
//~ cout << endl; 

		dL = std::max(dL, 0.);
		dxdt[Lind()] = dL; 
    }
    
    state_type generate_initial_conditions( int m){
		state_type x( (int)pow(m,2) + 1, 0. ); 
		int k = 0; 
		for (int i =0; i < m; i++){
			for (int j =0 ; j < m; j++){
				if (i == j){
					x[k] = p2logodds(1.); 
				} else{
					x[k] = p2logodds(0.);
				}
				k++; 
			}
		}
		return x; 
	}

	arma::mat Q_from_state(  state_type xfin ){
		arma::mat Q = arma::zeros(m,m); 
		int k =0 ;
		for (int i = 0; i < m; i++){
			for (int j = 0; j < m; j++){
				Q.at(i,j) = logodds2p( xfin[k] ); 
				k++; 
			}
		}
		for (int i = 0; i < m; i++){
			Q.row(i) = Q.row(i) / sum(Q.row(i)); 
		}
		return Q;
	}

	arma::vec A_from_state( state_type x)
	{
		arma::vec A =  Q_from_state( x ).t() * A0 ; 
		A = sumA * A / sum(A) ;  
		return A; 
	}

	double L_from_state( state_type xfin){
		return (double)xfin[xfin.size()-1]; 
	}

    
private:
	double dbl_clamp( double x, double lb, double ub){
		if (x < lb ) return lb; 
		if (x > ub ) return ub; 
		return x;
	}
	
	double p2logodds( double p ){
		double pp = dbl_clamp( p, MIN_P, MAX_P); 
		return log( pp / (1.-pp));
	}
	
	double logodds2p( double u){
		//~ return 1. / (1. + exp(-u) );
		return dbl_clamp( 1. / (1. + exp(-u) ), MIN_P, MAX_P);
	}
	
	double Q( const state_type &x, int k, int l) {
		return logodds2p( x[l*m+k] ); 
	}
	double L( const state_type &x ) {
		return x[ (int)pow(m,2) ] ;
	}
	int Qind( int k, int l ){
		return l*m + k ; 
	}
	int Lind(){
		return (int)pow(m,2)  ; 
	}
}; 




//[[Rcpp::export()]]
List solveQALboost1(arma::vec times, List Fs, List Gs, List Ys
 , double h0
 , double h1
 , double L0
 , arma::vec A0
){
	double treeT = std::abs( times(0) - times(times.size()-1)); 
	int m = A0.size();
	arma::mat Q0 = arma::diagmat( arma::ones( m )) ;
	double hres = times.size();
	
	DQAL3 dqal(Fs, Gs, Ys, m, hres, treeT, A0 ); 
	state_type x = dqal.generate_initial_conditions( A0.size() ); 
	
	//~ size_t steps = boost::numeric::odeint::integrate( dqal ,  x , h0 , h1 
			  //~ , (h1-h0)/10. );  
	
	boost::numeric::odeint::runge_kutta4< state_type > stepper;
	boost::numeric::odeint::integrate_const( stepper , dqal , x , h0, h1 ,  (h1-h0)/10. );
	
	arma::mat Q = dqal.Q_from_state( x); 
	arma::vec A = dqal.A_from_state( x); 
	double L = dqal.L_from_state(x); 
	
//~ cout << h0 << endl;
//~ cout << h1 << endl;
//~ cout << Q << endl;
//~ cout << A0 << endl;
//~ cout << A << endl;
//~ cout << L << endl;
//~ cout << endl;
//~ cout << endl;

	List o; 
	o["Q"] = Q;
	o["A"] = A;
	o["L"] = L + L0;
	return o; 
}
