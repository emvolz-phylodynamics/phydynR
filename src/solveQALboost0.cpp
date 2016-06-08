// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>

#include <RcppArmadillo.h>



using namespace arma;
using namespace Rcpp; 
using namespace std; 


typedef std::vector<double> state_type; 


// note: this is a simple copy/paste from colik.cpp for debugging. This should also match contents of rcolgem package;
class DQAL2{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
public:
	DQAL2( List Fs, List  Gs, List Ys, int m, double hres, double treeT ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT) {};
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
    {
		// time index
		//~ int i =  (int)min( 1+(int)( hres * (*t) / treeT ), hres);
		//~ NOTE hres = length(times)
		//~ NOTE treeT = max(times)
		int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); 
		mat F = as<mat>(Fs[i]); 
		mat G = as<mat>(Gs[i]); 
		vec Y = as<vec>(Ys[i]); 
		
		int k,l,z,w;
		
		double a[m]; //normalized nlft 
		//~ double sumA = 0.; 
		//~ for (k = 0; k < m; k++) sumA += A(k);
		double r = 1. ; // Atotal / sumA; // TODO
		for (k = 0; k < m; k++) { 
			//dA(k) = 0.;
			dxdt[Aind(k)] = 0.; 
			if (Y(k) > 0) {
				//~ a[k] = r *  std::min(1., A(x, k)/ Y(k) );
				a[k] = r *  A(x, k)/ Y(k) ;
			} else{
				a[k] = 1.; //
			} 
		}
		
		//dA
		for (k = 0; k < m; k++){
			for (l = 0; l < m; l++){
				if (k==l){
					//dA(k) -= a[l] * (F(i,l,k)) * a[k];
					dxdt[ Aind(k) ] -= a[l] * F(l,k) * a[k];
				} else{
					dxdt[ Aind(k) ]  += ( std::max(0., (1 - a[k])) * F(k,l) + G(k,l)) * a[l] ;
					dxdt[ Aind(k) ]  -= (F(l,k) + G(l,k)) * a[k];
				}
			}
		}
		//dQ
		for (z = 0; z < m; z++){ // col of Q
			for (k = 0; k < m; k++){ //row of Q
				dxdt[ Qind(k, z ) ] = 0. ; 
				for (l = 0. ; l < m; l++){
					if (k!=l){
						if ( Q(x, l,z) > 0)
						{
							dxdt[ Qind(k, z ) ] += (F(k,l) + G(k,l)) *  Q(x, l,z)/  std::max(Q(x,l,z), Y(l));
						}
						if (Q(x, k,z) > 0)
						{
							dxdt[ Qind(k, z ) ] -= (F(l,k) + G(l,k)) *  Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
						}
					}
					// coalescent:
					//~ dQ(k,z) -= (F(i,k,l)+F(i,l,k)) * a[l] * Q(k,z)/Y(i,k);
					if (Q(x, k,z) > 0){
						dxdt[ Qind(k, z ) ] -= F(k,l) * a[l] * Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
					}
				}
			}
		}
		//dL
		double dL = 0.;
		double Ydenom;
		for (k = 0; k < m; k++){
			for (l =0 ; l < m; l++){
				if (k==l){
					if (Y(k) > 1 && A(x,k) > 1){
						Ydenom = std::max( (Y(k)-1.),(r*A(x,k)-1.) );//
						if (Ydenom > 0)
						{
							dL += std::max(0.,std::min(1.,a[k])) * (  (r*A(x, k)-1)/Ydenom) * F(k,l);
						}
					} else{
						dL += std::max(0.,std::min(1.,a[k])) * std::max(0.,std::min(1.,a[l])) * F(k,l);
					}
				} else{
					dL += std::max(0.,std::min(1.,a[k])) * std::max(0.,std::min(1.,a[l])) * F(k,l);
				}
			}
		}
		dL = std::max(dL, 0.);
		dxdt[Lind()] = dL; 
    }
    
    state_type generate_initial_conditions( vec A){
		int m = A.size() ;
		state_type x( (int)pow(m,2) + m + 1, 0. ); 
		int k = 0; 
		for (int i =0; i < m; i++){
			for (int j =0 ; j < m; j++){
				if (i == j){
					x[k] =1.; 
				} 
				k++; 
			}
		}
		for (int i = 0; i < m; i++){
			x[k] = A.at(i) ; 
			k++;
		}
		return x; 
	}

	mat Q_from_state(  state_type xfin ){
		mat Q = zeros(m,m); 
		int k =0 ;
		for (int i = 0; i < m; i++){
			for (int j = 0; j < m; j++){
				Q.at(i,j) = xfin[k]; 
				k++; 
			}
		}
		for (int i = 0; i < m; i++){
			Q.row(i) = Q.row(i) / sum(Q.row(i)); 
		}
		return Q;
		//~ Q = Q / arma::sum(Q, 1) ; 
	}

	vec  A_from_state( state_type xfin){
		vec A = zeros(m); 
		int k = 0; 
		for ( int i = (int)pow(m,2); i < ((int)pow(m,2)+m); i++){
			A(k) = xfin[i]; 
			k++; 
		}
		return A; 
	}

	double L_from_state( state_type xfin){
		return (double)xfin[xfin.size()-1]; 
	}

    
private:
	double Q( const state_type &x, int k, int l) {
		return x[l*m+k]; 
	}
	double A( const state_type &x, int k) {
		return x[(int)pow(m,2) + k]; 
	}
	double L( const state_type &x ) {
		return x[ (int)pow(m,2) + m] ;
	}
	
	int Qind( int k, int l ){
		return l*m + k ; 
	}
	int Aind( int k ){
		return (int)pow(m,2) + k; 
	}
	int Lind(){
		return (int)pow(m,2) + m ; 
	}
}; 



//[[Rcpp::export()]]
List solveQALboost0(vec times, List Fs, List Gs, List Ys
 , double h0
 , double h1
 , double L0
 , vec A0
 , double treeT
){
	int m = A0.size();
	mat Q0 = diagmat( ones( m )) ;
	double hres = times.size();//abs(times(1) - times(0)); 
	DQAL2 dqal(Fs, Gs, Ys, m, hres, treeT ); 
	state_type x = dqal.generate_initial_conditions( A0 ); 
//~ cout << hres << endl;
//~ Rf_PrintValue( wrap(x)); 
//~ cout << A0 << endl;
//~ cout << Q0 << endl;
//~ cout << endl;
	size_t steps = boost::numeric::odeint::integrate( dqal ,  x , h0 , h1 
			  , (h1-h0)/10. );  
	mat Q = dqal.Q_from_state( x); 
	vec A = dqal.A_from_state( x); 
	double L = dqal.L_from_state(x); 
//~ cout << Q << endl;
//~ cout << A << endl;
//~ cout << L << endl;
	List o; 
	o["Q"] = Q;
	o["A"] = A;
	o["L"] = L;
	return o; 
}
