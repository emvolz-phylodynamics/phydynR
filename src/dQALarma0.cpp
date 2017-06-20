/*
 * derivative function to be used with deSolve
 */ 
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>



using namespace arma;
using namespace Rcpp; 
using namespace std; 


const double MIN_Y = 1. ;

mat Q_from_state( const vec xfin, int m ){
	// TODO faster version using fill
	mat Q = zeros(m,m); 
	int k =0 ;
	for (int i = 0; i < m; i++){
		for (int j = 0; j < m; j++){
			Q(i,j) = xfin[k]; 
			k++; 
		}
	}
	for (int i = 0; i < m; i++){
		Q.row(i) = Q.row(i) / sum(Q.row(i)); 
	}
	return Q;
}

int Qind( int k, int l, int m ){
	return l*m + k ; 
}
int Lind(int m){
	return (int)pow(m,2)  ; 
}

// replaced all Q(k,z) with Q(x,k,z,m)
double Q( const vec &x, int k, int l, int m) {
     return x[l*m+k]; 
}

//[[Rcpp::export()]]
vec dQL( vec x , mat F , mat G, vec Y, vec A0){
	int m = Y.size();
	Y = clamp(Y, MIN_Y, INFINITY);
	vec dx = zeros<vec>(m*m+1); 
	// mat Q = Q_from_state(x, m );  -- use Q(..) instead
	//vec A = Q.t() * A0;
	
        vec A =  Q_from_state( x,m ).t() * A0 ; 
	A = sum(A0) * A / sum(A);  // as done in A_from_state
	
	vec a = A / Y ;
	
	int Lind = x.size();
	int z, k, l ;
	//dQ
	for (z = 0; z < m; z++){ // col of Q
	   for (k = 0; k < m; k++){ //row of Q
	      for (l = 0. ; l < m; l++){
		if (k!=l){
		  if ( Q(x,l,z,m) > 0) {
		    dx(Qind(k, z,m)) += (F(k,l) + G(k,l)) *  Q(x,l,z,m)/ std::max(Q(x,l,z,m), Y(l));
		  }
		  if (Q(x,k,z,m) > 0) {
		    dx( Qind(k, z, m )) -= (F(l,k) + G(l,k)) *  Q(x,k,z,m)/  std::max(Q(x,k,z,m), Y(k));
		  }
		}
	       // coalescent: //TODO this should use a slightly different Q
		if (Q(x,k,z,m) > 0){
		  dx( Qind(k, z,m ) ) -= F(k,l) * a(l) * Q(x,k,z,m)/  std::max(Q(x,k,z,m), Y(k));
	       }
	     }
	   }
	}
	//dL
	double dL = 0.;
	for (k = 0; k < m; k++){
	   for (l =0 ; l < m; l++){				
	      if (k == l && A(k) >= 1. ){
		dL += (A(k) / Y(k)) * ((A(k)-1.) / Y(k)) * F(k,l) ; 
	      } else {
		dL += a(k) * a(l) * F(k,l);
	      }
	    }
	}
	dx(dx.size()-1) = dL;
	
	return dx; 
}
