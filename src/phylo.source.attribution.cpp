/* lineages through time and hazard of co
 */


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>
#include <RcppArmadillo.h>
#include "colik.hpp" 


using namespace arma;
using namespace Rcpp; 
using namespace std; 
using namespace boost::numeric::odeint; 

//~ static const int SAMPLE  = 0; 
//~ static const int CO = 1; 
static const double  MIN_Y = 1.; 


// [[Rcpp::export]]
NumericMatrix updateWCpp(NumericMatrix W,  NumericVector psi_a
  ,  IntegerVector utips
  ,  IntegerVector vtips
  ,  IntegerVector utips_Wcoords
  ,  IntegerVector vtips_Wcoords)
{
	//~ NumericMatrix WW  = clone( W);
	int ut,vt,utW,vtW;
	for (int iu = 0; iu < utips.size(); iu++){
		for (int iv = 0; iv < vtips.size(); iv++){
			ut = utips(iu)-1;
			vt = vtips(iv)-1;
			utW = utips_Wcoords(iu)-1;
			vtW = vtips_Wcoords(iv)-1;  
			//~ WW(utW,vtW) = WW(utW,vtW) + psi_a(ut) * psi_a(vt) / 2.; 
			//~ WW(vtW,utW) = WW(utW,vtW); 
			W(utW,vtW) = W(utW,vtW) + psi_a(ut) * psi_a(vt) / 2.; 
			W(vtW,utW) = W(utW,vtW); 
		}
	}
	return W; 
}


///////////////////////////



typedef std::vector<double> ode_state_type; 

class DPikRhoPsi{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
	vec A0; 
	double sumA; 
	int nextant; 
	int nextant_tip; 
public:
	DPikRhoPsi( const List Fs, const List  Gs, const List Ys, const int m, const double hres, const double treeT, const int nextant) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT),nextant(nextant),nextant_tip(nextant_tip) { };
	
	void operator() ( const ode_state_type &x , ode_state_type &dxdt , double t)//, const double /* t */ )
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
		
		mat R_rho = G; 
		R_rho.each_row() /= Y.t();
		R_rho.diag().zeros(); 
		R_rho.diag() = R.diag(); // NOTE prob not conserved, sum should be psi //-sum( R, 0).t(); 
		
		mat Pik = x2Pik( x ); 
		mat dPik = zeros<mat>( m, nextant); 
		
		mat rho_ik = x2rho( x ); 
		mat drho_ik = zeros<mat>( m, nextant_tip); 
		
		vec A = sum(Pik, 1 );
		vec Ampik = zeros(m); 
		vec Ampik_Y = zeros(m); 
		// d P
		for (z = 0; z < nextant; z++){
			dPik.col(z) = R * Pik.col(z) ; 
			
			Ampik = clamp( A - Pik.col(z), 0., INFINITY);

			Ampik_Y = clamp( A / Y , 0., 1e6 );
			for (k = 0; k < m; k++){
				for( l = 0; l < m; l++){
					dPik(k,z) -= (Pik(k,z)/Y(k)) * F(k,l) * Ampik_Y(l) ; 
				}
			}
		}
		
		// d rho
		for (z = 0; z < nextant_tip; z++){
			drho_ik.col(z) = R_rho * rho_ik.col(z) ; 
		}
		
		//pass
		w = 0; 
		for (int z = 0; z < nextant; z++){
			for (int k = 0; k < m; k++){
				dxdt[w] = dPik( k, z); 
				w++;
			}
		}
		
		//~ w = 0; 
		for (int z = 0; z < nextant_tip; z++){
			for (int k = 0; k < m; k++){
				dxdt[w] = drho_ik( k, z); 
				w++;
			}
		}
	}
	
	mat x2Pik( ode_state_type x ){
		return mat( &x[0], m, nextant ); 
	}
	
	mat x2rho( ode_state_type x ){
		return mat( &x[m*nextant], m, nextant_tip ); 
	}
	
	void set_nextant(int ne ){ 
		nextant = ne; 
	}
	
	void set_nextant_tip(int ne ){ 
		nextant_tip = ne; 
	}
	
	ode_state_type generate_initial_conditions( mat pik0, mat rho0, int m, int nextant, int nextant_tip ){
		ode_state_type x( m * nextant + m * nextant_tip , 0.); 
		int w = 0; 
		for (int z =0; z < nextant; z++){
			for (int k = 0; k < m; k++){
				x[w] = abs(pik0(k,z));
				w++; 
			}
		}
		for (int z = 0; z < nextant_tip; z++){
			for ( int k = 0; k < m; k++){
				x[w] = abs(rho0(k,z)); 
				w++;
			}
		}
		return x;
	}
};


//[[Rcpp::export()]]
List sourceAttribMultiDemeCpp2( const NumericVector heights, const List Fs, const List Gs, const List Ys
  , const IntegerVector eventIndicator // sample or co
  , const IntegerVector eventIndicatorNode // node involved at each event
  , const NumericVector eventHeights
  , const mat sampleStates  
  , const IntegerMatrix daughters // daughters of each node
  , const int n
  , const int Nnode
  , const int m
  , double AgtYboundaryCondition
  , const double maxHeight // terminate computation at this height
  , const int step_size_res 
){
	double loglik = 0.;
	mat P(m, n + Nnode, fill::zeros); 
	mat rho(m, n , fill::zeros); 
	int u,v,w,z,a, k, l; 
	int samplesAdded = 0;
	vec psi = zeros(m) ; // corresponding to each state over each interval
	vec psi_time = zeros(n); // for each tip, varies over time 
	
	vec A_Y;
	vec Y ;
	vec A = zeros(m);
	mat F(m,m, fill::zeros);
	mat G(m,m, fill::zeros); 
	vec puY, pvY, pa;
	vec m_rs_R; 
	//~ std::vector<bool> extant(n + Nnode, false); 
	uvec extant = zeros<uvec>(n+Nnode);
	uvec extant_tip = zeros<uvec>(n); 
	uvec extant_indices;
	uvec extant_tip_indices;
	int nextant = 0;  
	int nextant_tip = 0;
	
	umat tipsDescendedFrom(n+Nnode, n, fill::zeros);
	for (u = 0; u < n; u++){
		tipsDescendedFrom.at(u,u) = 1; 
	} 
	// container for output 
	//~ mat W = zeros(n,n); 
	std::vector<int> donorW; 
	donorW.reserve( (int)(n * std::sqrt((double)n) ) ); 
	std::vector<int> recipW; 
	recipW.reserve( (int)(n * std::sqrt((double)n) ) ); 
	std::vector<double> W; 
	W.reserve( (int)(n * std::sqrt((double)n) ) ); 
	
	// instantiate solver 
	double hres =  heights.size() ;
	double treeT = heights[heights.size()-1]; // note increasing order   
	
	double psi_w, psi_z;   
	
	ode_state_type x;
	DPikRhoPsi dpikrho(Fs, Gs, Ys, m, hres, treeT, nextant ); 
	runge_kutta_cash_karp54<ode_state_type> stepper;
	double dt = std::abs( heights(1) - heights(0));
		
	// iterate events 
	double nextEventHeight = eventHeights(0); 
	int ievent = 0; 
	
	double sterm0, sterm1; 
	double h, h0, h1, hstar, dh; 
	int ih; // index of nextEventHeight
	size_t integrate_steps;
	
	//
	uvec row_indices = zeros<uvec>(m); 
	for (int k = 0;k < m; k++) row_indices(k) = k; 
	
	h = 0.; 
	while( nextEventHeight != INFINITY && nextEventHeight < maxHeight ){
		extant_indices = find(extant);
		extant_tip_indices = find(extant_tip);
		dpikrho.set_nextant( nextant ); 
		dpikrho.set_nextant_tip( nextant_tip ); 
		
		if (nextEventHeight > h )
		{
			A = sum(P, 1);
			A = arma::normalise(A,1.)  * ((double)nextant);
			x = dpikrho.generate_initial_conditions( P.submat(row_indices, extant_indices)
			  , rho.submat(row_indices, extant_tip_indices)
			  , m, nextant, nextant_tip); 
			//
			if (step_size_res < 1){
				integrate_steps = integrate( dpikrho ,  x , h , nextEventHeight
						  , (nextEventHeight - h)/100. );  
			} else{ 
				integrate_const( stepper, dpikrho, x, h, nextEventHeight
				 ,  std::min( dt/step_size_res, (nextEventHeight-h)/step_size_res) );
			}

			P.submat(row_indices, extant_indices) = dpikrho.x2Pik( x );
			rho.submat(row_indices, extant_tip_indices) = dpikrho.x2rho( x );
			
		} 
		
		if (eventIndicator(ievent)==SAMPLE)
		{
			u = eventIndicatorNode(ievent); 
			nextant++; 
			nextant_tip++; 
			//~ P.col(u-1) = arma::normalise(sortedSampleStates.col(samplesAdded) ,1.); //TODO
			P.col(u-1) = arma::normalise(sampleStates.col(u-1) ,1.);
			rho.col(u-1) = P.col(u-1); 
			psi_time.at(u-1) = 1.; 
			extant(u-1) = 1; 
			extant_tip(u-1) = 1; 
			samplesAdded++; 
		} else {
			// coalescent
			ih = (int)std::min( hres * nextEventHeight / treeT, (double)(heights.size()-1)); 
			F = as<mat>(Fs[ih]);
			G = as<mat>(Gs[ih]);
			Y = as<vec>(Ys[ih]);
			Y = arma::clamp(Y, 1e-6, INFINITY );
			
			a = eventIndicatorNode(ievent);
			u = daughters(a -1 , 0); 
			v = daughters( a - 1, 1 ); 
			puY = arma::normalise(  arma::min(Y,P.col(u- 1 )) ,1.) / arma::clamp( Y, MIN_Y, INFINITY ) ; 
			pvY = arma::normalise(  arma::min(Y,P.col(v- 1 )) ,1.) / arma::clamp( Y, MIN_Y, INFINITY ) ; 
			
			// state of ancestor 
			pa =  arma::normalise( (F * puY) % pvY + (F * pvY) % puY ,1.) ; 
			//~ pa = pa / sum(pa ) ; 
			P.col(a - 1 ) = pa; 
			
			//bookkeeping
			nextant--; 
			extant(u-1) = 0;
			extant(v-1) = 0;
			extant(a-1) = 1;
			P.col(u-1) = zeros<colvec>(P.n_rows); 
			P.col(v-1) = zeros<colvec>(P.n_rows);
			
			// sa stuff ; upate psi , rho 
			vec rho_w__Y = zeros(m); 
			vec rho_z__Y = zeros(m); 
			double pwz0,pzw0,pwz,pzw;
			// update W(iw, iz)
			for (int iw = 0; iw < n; iw++){
				if (tipsDescendedFrom.at(u-1, iw)==1){
					rho_w__Y = arma::normalise(  arma::min(Y,rho.col(iw)) ,1.) / arma::clamp(Y, MIN_Y, INFINITY ) ; 
					for (int iz = iw+1; iz < n; iz++){
						if (tipsDescendedFrom.at(v-1, iz)==1){
							rho_z__Y = arma::normalise(  arma::min(Y,rho.col(iz)) ,1.) / arma::clamp(Y, MIN_Y, INFINITY ) ; 
							pwz0 = sum( rho_w__Y % (F * rho_z__Y) );
							pzw0 = sum( rho_z__Y % (F * rho_w__Y ));
							
							psi_w = sum( rho.col(iw)); 
							psi_z = sum(rho.col(iz)); 
							if ( (pwz0 + pzw0) > 0 ){
								pwz = psi_w * psi_z * pwz0 / ( pwz0 + pzw0); 
								pzw = psi_w * psi_z * pzw0 / ( pwz0 + pzw0); 
								
							}  else{
								//~ cout << " (pwz0 + pzw0) = 0 " << endl; 
								//~ cout << psi_time.at(iw) << " " <<  psi_time.at(iz)<< endl; 
								//~ cout << rho.col(iw) << endl; 
								//~ cout << rho.col(iz) << endl;
								pwz = 0.;
								pzw = 0.; 
							}
							
							donorW.push_back( iw + 1);
							recipW.push_back( iz+1 );
							W.push_back( pwz) ;
							
							donorW.push_back( iz + 1 );
							recipW.push_back( iw+ 1);
							W.push_back( pzw );
							
							//~ W.at( iw, iz) = pwz; 
							//~ W.at(iz, iw) = pzw; 
						}
					}
				}
			}
			tipsDescendedFrom.row(a-1) = tipsDescendedFrom.row(u-1) + tipsDescendedFrom.row(v-1); 
			//update rho and psi
			for (int iw = 0; iw < n; iw++){// w & z tips
				if (tipsDescendedFrom.at(u-1, iw)==1){
					//tip descended from a and u
					rho_w__Y = arma::normalise( arma::clamp(  arma::min(Y,rho.col(iw))  / Y, 0., 1. ) ,1.) ; 
					//update psi(iw)
					double puv = sum( rho_w__Y % (F * pvY) );
					double pvu = sum( pvY % (F * rho_w__Y ));
					puv = puv / (puv + pvu ) ;
					psi_time.at(iw) = sum( rho.col(iw)); 
					psi_time.at(iw) *= puv ;
					
					//update rho(iw)
					rho.col(iw) =  (rho_w__Y % (F * pvY)) ;
					rho.col(iw) =  psi_time.at(iw) * rho.col(iw) / sum(rho.col(iw));
				} else if (tipsDescendedFrom.at(v-1, iw)==1){
					//tip descended from a and v
					rho_w__Y = arma::normalise(arma::clamp(  arma::min(Y,rho.col(iw))  / Y, 0., 1. ) ,1.) ; 
					//update psi(iw)
					double pvu = sum( rho_w__Y % (F * puY) );
					double puv = sum( puY % (F * rho_w__Y ));
					pvu = pvu / (puv + pvu ) ;
					psi_time.at(iw) = sum( rho.col(iw)); 
					psi_time.at(iw) *= pvu ;
					//update rho(iw)
					rho.col(iw) =  (rho_w__Y % (F * pvY)) ;
					rho.col(iw) =  psi_time.at(iw) * rho.col(iw) / sum(rho.col(iw));
				}
			}
		}
		h = nextEventHeight; 
		ievent++;
		if (ievent<eventHeights.size()){
			nextEventHeight = eventHeights(ievent);
		} else{
			nextEventHeight = INFINITY; 
		}
		if ( h > maxHeight ) {
			break; 
		}
	}
	List wlist; 
	wlist["donor"] = wrap(donorW);
	wlist["recip"] = wrap(recipW);
	wlist["infectorProbability"] = wrap(W);
	return wlist; 
}


