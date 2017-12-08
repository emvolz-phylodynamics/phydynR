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

typedef std::vector<double> state_type; 



class DQfwd{
	List Fs, Gs, Ys, deaths; 
	int m; 
	double hres; 
	double treeT; 
	arma::vec A0; 
	double sumA; 
public:
	DQfwd( List Fs, List  Gs, List Ys, List deaths, int m, double hres, double treeT ) : Fs(Fs),Gs(Gs),Ys(Ys),deaths(deaths),m(m),hres(hres),treeT(treeT){ };
	
	
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
    {
		// time index-- NOTE time is on forward axis
		// FGYdeaths is on reverse axis
		//~ NOTE hres = length(times)
		//~ NOTE treeT = max(times)
		int i =  (int)std::max(0., std::min( hres * (treeT-t) / treeT , (double)(hres-1.))); 
		arma::mat F = as<arma::mat>(Fs[i]); 
		arma::mat G = as<arma::mat>(Gs[i]); 
		arma::vec Y = arma::clamp(as<arma::vec>(Ys[i]), MIN_Y, INFINITY ); 
		arma::vec death = as<arma::vec>(deaths[i]); 
		
		int k,l,z,w;
		
		//dQ
		for (z = 0; z < m; z++){ // col of Q
			for (k = 0; k < m; k++){ //row of Q
				dxdt[ Qind(k, z ) ] = 0. ; 
				for (l = 0. ; l < m; l++){
					if (k!=l){
						if ( Q(x, l,z) > 0)
						{
							dxdt[ Qind(k, z ) ] += (F(l,k)/2. + G(l,k)) *  Q(x,l,z)/  std::max(Q(x,l,z), Y(l));
						}
						if (Q(x, k,z) > 0)
						{
							dxdt[ Qind(k, z ) ] -= (F(k,l)/2. + G(k,l)) *  Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
						}
					}
				}
				// death:
				if (Q(x, k,z) > 0){
					dxdt[ Qind(k, z ) ] -= Q(x, k,z) * death(k) / std::max(Q(x, k,z), Y(k));
				}
			}
		}
    }
    
    state_type generate_initial_conditions( int m){
		state_type x( (int)pow(m,2) , 0. ); 
		int k = 0; 
		for (int i =0; i < m; i++){
			for (int j =0 ; j < m; j++){
				if (i == j){
					x[k] =1.; 
				} 
				k++; 
			}
		}
		return x; 
	}
	
	// note returns Q in col order
	arma::mat Q_from_state(  state_type xfin ){
		arma::mat QQ = arma::zeros(m,m); 
		int k =0 ;
		for (int j = 0; j < m; j++){
			for (int i = 0; i < m; i++){
				QQ.at(i,j) = xfin[k]; 
				k++; 
			}
		}
		for (int i = 0; i < m; i++){
			QQ.col(i) = QQ.col(i) / sum(QQ.col(i)); 
		}
		return QQ;
	}
private:
	double Q( const state_type &x, int k, int l) {
		return x[l*m+k]; 
	}
	int Qind( int k, int l ){
		return l*m + k ; 
	}
}; 




//[[Rcpp::export()]]
arma::mat solveQfwd0(arma::vec times, List Fs, List Gs, List Ys, List deaths
 , int m
 , double h1
 , double h0
){
	double treeT = std::abs( times(0) - times(times.size()-1)); 
	arma::mat Q0 = arma::diagmat( arma::ones( m )) ;
	double hres = times.size();//abs(times(1) - times(0)); 
	
	DQfwd dqfwd  (Fs, Gs, Ys, deaths, m, hres, treeT ); 
	state_type x = dqfwd.generate_initial_conditions( m ); 
	double t0 = treeT - h1; 
	double t1 = treeT - h0; 
	size_t steps = boost::numeric::odeint::integrate( dqfwd ,  x , t0 , t1 
			  , (h1-h0)/100. );  
	
	arma::mat Q = dqfwd.Q_from_state( x); 
	return Q; 
}
