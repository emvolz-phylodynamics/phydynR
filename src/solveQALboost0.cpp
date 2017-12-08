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
		arma::mat F = as<arma::mat>(Fs[i]); 
		arma::mat G = as<arma::mat>(Gs[i]); 
		arma::vec Y = as<arma::vec>(Ys[i]); 
		
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
							dxdt[ Qind(k, z ) ] -= (F(l,k) + G(l,k)) *  Q(x, k,z)/  std::max(Q(x,k,z), Y(k));
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
    
    state_type generate_initial_conditions( arma::vec A){
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

	arma::mat Q_from_state(  state_type xfin ){
		arma::mat Q = arma::zeros(m,m); 
		int k =0 ;
		for (int i = 0; i < m; i++){
			for (int j = 0; j < m; j++){
				Q.at(i,j) = xfin[k]; 
				k++; 
			}
		}
		for (int i = 0; i < m; i++){
			Q.row(i) = Q.row(i) / arma::sum(Q.row(i)); 
		}
		return Q;
		//~ Q = Q / arma::sum(Q, 1) ; 
	}

	arma::vec  A_from_state( state_type xfin){
		arma::vec A = arma::zeros(m); 
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


// this version does not solve extra eqns for A; also conserves sum A_k
class DQAL2_2{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
	arma::vec A0; 
	double sumA; 
public:
	DQAL2_2( List Fs, List  Gs, List Ys, int m, double hres, double treeT, arma::vec _xA0 ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT),A0(_xA0) {
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
		arma::vec Y = clamp(as<arma::vec>(Ys[i]), MIN_Y, INFINITY ); 
		
		int k,l,z,w;
		
		arma::vec A = A_from_state( x ) ;
		//~ vec a = clamp( A / Y , 0., 1. ); 
		arma::vec a = A / Y ; 
				
//~ cout << i << endl;
//~ cout << Y << endl;
//~ cout << a << endl;
//~ cout  << endl;
//~ cout  << endl;

		
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
					if (Q(x, k,z) > 0){
						dxdt[ Qind(k, z ) ] -= F(k,l) * a(l) * Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
					}
				}
			}
		}
		//dL
		double dL = 0.;
		double Ydenom;
		for (k = 0; k < m; k++){
			for (l =0 ; l < m; l++){				
				if (k == l && A(k) >= 1. ){
					//~ dL += (A(k) / Y(k)) * ((A(k)-1.) / Y(k)) * F(k,l); 
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
					x[k] =1.; 
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
				Q.at(i,j) = xfin[k]; 
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
	double Q( const state_type &x, int k, int l) {
		return x[l*m+k]; 
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
List solveQALboost0(arma::vec times, List Fs, List Gs, List Ys
 , double h0
 , double h1
 , double L0
 , arma::vec A0
 //, double treeT
){
	double treeT = std::abs( times(0) - times(times.size()-1)); 
	int m = A0.size();
	arma::mat Q0 = arma::diagmat( arma::ones( m )) ;
	double hres = times.size();//abs(times(1) - times(0)); 
	
	//~ DQAL2 dqal(Fs, Gs, Ys, m, hres, treeT ); 
	//~ state_type x = dqal.generate_initial_conditions( A0 ); 
	DQAL2_2 dqal(Fs, Gs, Ys, m, hres, treeT, A0 ); 
	state_type x = dqal.generate_initial_conditions( A0.size() ); 
	
	size_t steps = boost::numeric::odeint::integrate( dqal ,  x , h0 , h1 
			  , (h1-h0)/100. );  
	
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
