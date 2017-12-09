#ifndef __COLIK_HPP_INCLUDED__
#define __COLIK_HPP_INCLUDED__



#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//~ TODO should be dfined here 
static const int SAMPLE  = 0; 
static const int CO = 1; 
typedef std::vector<double> ode_state_type; 

mat finite_size_correction2(const vec& p_a, const vec& A, const std::vector<bool> extant, mat P);

class DQAL{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
public:
	DQAL( List Fs, List  Gs, List Ys, int m, double hres, double treeT );
	void operator() ( const ode_state_type &x , ode_state_type &dxdt , double t);
private:
	double Q( const ode_state_type &x, int k, int l);
	double A( const ode_state_type &x, int k);
	double L( const ode_state_type &x );
	
	int Qind( int k, int l );
	int Aind( int k );
	int Lind();
}; 


ode_state_type generate_initial_conditions( vec A);

void Q_from_state( mat &Q, ode_state_type xfin );

void A_from_state( vec &A, ode_state_type xfin);

double L_from_state( ode_state_type xfin);

double colik2cpp(const NumericVector heights, const List Fs, const List Gs, const List Ys
  , const IntegerVector eventIndicator // sample or co
  , const IntegerVector eventIndicatorNode // node involved at each event
  , const NumericVector eventHeights
  , const mat sortedSampleStates
  , const IntegerMatrix daughters // daughters of each node
  , const int n
  , const int Nnode
  , const int m
  , bool AgtYboundaryCondition );


#endif
