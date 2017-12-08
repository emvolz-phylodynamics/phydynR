// [[Rcpp::depends(RcppArmadillo, BH)]]

#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/euler.hpp>

#include <RcppArmadillo.h>



using namespace arma;
using namespace Rcpp; 

static const int SAMPLE  = 0; 
static const int CO = 1; 

typedef std::vector<double> ode_state_type; 

mat finite_size_correction2(const vec& p_a, const vec& A, const std::vector<bool> extant, mat P)
{
	// NOTE mstates m X n 
	int u; 
	vec rho; 
	vec rterm; 
	//~ vec lterm; 
	double lterm; 
	vec p_u; 
	vec Amin_pu; 
	//~ for (int iu = 0; iu < extantLines.size(); iu++){
	for (int iu = 0; iu < extant.size(); iu++){
		if (extant.at(iu)){
			u = iu + 1; 
			p_u = P.col(u-1); 
			Amin_pu = clamp(( A - p_u), 1., INFINITY ); 
			rterm = p_a / Amin_pu ; 
			rho = A / Amin_pu; 
			lterm = dot( rho, p_a); //
			p_u = p_u % (lterm - rterm) ;
			p_u = p_u / sum(p_u ) ; 
			P.col(u - 1) = p_u; 
		}
	}
	return P; 
}

void finite_size_correction3(const int ia, const vec& A, const std::vector<bool> extant, mat& P)
{
	// NOTE mstates m X n 
	int u; 
	vec rho; 
	vec rterm; 
	//~ vec lterm; 
	double lterm; 
	vec p_u; 
	vec Amin_pu; 
	vec p_a = P.col( ia ); 
	for (int iu = 0; iu < extant.size(); iu++){
		if (extant.at(iu) && iu!=ia){
			u = iu + 1; 
			p_u = P.col(u-1); 
			Amin_pu = clamp(( A - p_u), 1., INFINITY ); 
			rterm = p_a / Amin_pu ; 
			rho = A / Amin_pu; 
			lterm = dot( rho, p_a); //
			p_u = p_u % clamp((lterm - rterm), 0., INFINITY ) ; // l > r
			p_u = p_u / sum(p_u ) ; 
			P.col(u - 1) = p_u; 
		}
	}
}


class DQAL{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
public:
	DQAL( List Fs, List  Gs, List Ys, int m, double hres, double treeT ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT) {};
	void operator() ( const ode_state_type &x , ode_state_type &dxdt , double t)//, const double /* t */ )
    {
		// time index
		//~ int i =  (int)min( 1+(int)( hres * (*t) / treeT ), hres);
		//~ NOTE hres = length(times)
		//~ NOTE treeT = max(times)
		int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); //TODO redefine max min? will use arma versin? 
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
		// TODO shouldn't sum(A) be conserved over each interval?? then these are not the right eqns...
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
    
private:
	//~ double F(int ih, int k, int l){
		//~ return as<mat>(Fs[ih]).at(k,l);
	//~ }
	//~ double G(int ih, int k, int l){
		//~ return as<mat>(Gs[ih]).at(k,l);
	//~ }
	//~ double Y(int ih, int k){
		//~ return as<vec>(Ys[ih]).at(k);
	//~ }
	
	double Q( const ode_state_type &x, int k, int l) {
		return x[l*m+k]; 
	}
	double A( const ode_state_type &x, int k) {
		return x[(int)pow(m,2) + k]; 
	}
	double L( const ode_state_type &x ) {
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


ode_state_type generate_initial_conditions( vec A){
	int m = A.size() ;
	ode_state_type x( (int)pow(m,2) + m + 1, 0. ); 
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

void Q_from_state( mat &Q, ode_state_type xfin ){
	int m = Q.n_rows; 
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
	//~ Q = Q / arma::sum(Q, 1) ; 
}

void A_from_state( vec &A, ode_state_type xfin){
	int m = A.size(); 
	int k = 0; 
	for ( int i = (int)pow(m,2); i < ((int)pow(m,2)+m); i++){
		A.at(k) = xfin[i]; 
		k++; 
	}
}

double L_from_state( ode_state_type xfin){
	return (double)xfin[xfin.size()-1]; 
}




////////////////////////////////////////////////////////////////////////////////


// transition probabilities and survivor functions, conditioning on lineage not changing host
class DQpsiA{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
public:
	DQpsiA( List Fs, List  Gs, List Ys, int m, double hres, double treeT ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT) {};
	void operator() ( const ode_state_type &x , ode_state_type &dxdt , double t)//, const double /* t */ )
    {
		// time index
		//~ NOTE hres = length(times)
		//~ NOTE treeT = max(times)
		int i =  (int)std::max(0., std::min( hres * t / treeT , (double)(hres-1.))); 
		mat F = as<mat>(Fs[i]); 
		mat G = as<mat>(Gs[i]); 
		vec Y = as<vec>(Ys[i]); 
		
		int k,l,z,w;
		
		double a[m]; //normalized nlft 
		for (k = 0; k < m; k++) { 
			dxdt[Aind(k)] = 0.; 
			if (Y(k) > 0) {
				a[k] = A(x, k)/ Y(k) ;
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
		// NOTE this is filling role of "rho" in the paper 
		for (z = 0; z < m; z++){ // col of Q
			for (k = 0; k < m; k++){ //row of Q
				dxdt[ Qind(k, z ) ] = 0. ; 
				for (l = 0 ; l < m; l++){
					if (k!=l){
						if ( Q(x, l,z) > 0)
						{
							dxdt[ Qind(k, z ) ] += ( G(k,l)) *  Q(x, l,z)/  std::max(Q(x,l,z), Y(l));
						}
						if (Q(x, k,z) > 0)
						{
							dxdt[ Qind(k, z ) ] -= ( G(l,k)) *  Q(x, k,z)/  std::max(Q(x, k,z), Y(k));
						}
					}
				}
			}
		}
		// d psi  
		for (z = 0; z < m; z++){ // col of Q, elem of psi- line starts in deme z
			dxdt[ psi_ind( z ) ] = 0. ;
			for (k = 0; k < m; k++){ //row of Q - current deme
				for (l = 0. ; l < m; l++){ //source of transm
					// NOTE Q is filling the role of rho_{ik} in the paper 
					dxdt[psi_ind(z)] -= std::max(0., (1 - a[l])) * F(l,k) * std::max(0., psi(x, z)) * Q(x,k,z) / std::max(Q(x, k, z), Y(k)); 
				}
			}
		}
    }
private:
	double Q( const ode_state_type &x, int k, int l) {
		return x[l*m+k]; 
	}
	double psi( const ode_state_type &x, int k) {
		return x[(int)pow(m,2) +  k]; 
	}
	double A( const ode_state_type &x, int k) {
		return x[(int)pow(m,2) + m + k]; 
	}
	int Qind( int k, int l ){
		return l*m + k ; 
	}
	int psi_ind( int k ){
		return (int)pow(m,2) + k; 
	}
	int Aind( int k ){
		return (int)pow(m,2) + m + k; 
	}
};



ode_state_type generate_initial_conditions_dqpsia(vec A){
	int m = A.size() ;
	ode_state_type x( (int)pow(m,2) + m + m, 0. ); 
	int k = 0; 
	for (int i =0; i < m; i++){
		for (int j =0 ; j < m; j++){
			if (i == j){
				x[k] =1.; // Q or 'rho'
			} 
			k++; 
		}
	}
	for (int i = 0; i  < m ; i++){
		x[k] = 1.; //psi
		k++; 
	}
	for (int i = 0; i < m; i++){
		x[k] = A.at(i) ; //A
		k++;
	}
	return x; 
}


void Qrho_from_state( mat &Qrho, ode_state_type xfin ){
	int m = Qrho.n_rows; 
	int k =0 ;
	for (int i = 0; i < m; i++){
		for (int j = 0; j < m; j++){
			Qrho.at(i,j) = std::min(1.,std::max(0., xfin[k]));  
			k++; 
		}
	}
	for (int i = 0; i < m; i++){
		Qrho.row(i) = Qrho.row(i) / sum(Qrho.row(i)); 
	}
}


void psi_from_state( vec &psi, ode_state_type xfin ){
	int m = psi.size(); 
	int k = 0; 
	for ( int i = (int)pow(m,2); i < ((int)pow(m,2)+m); i++){
		psi.at(k) = std::min(1.,std::max(0., xfin[i])); 
		k++; 
	}
}


//[[Rcpp::export()]]
List sourceAttribMultiDemeCpp( const NumericVector heights, const List Fs, const List Gs, const List Ys
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
){
	double loglik = 0.;
	mat P(m, n + Nnode, fill::zeros); 
	mat rho(m, n , fill::zeros); 
	int nextant = 0; 
	int u,v,w,z,a, k, l; 
	int samplesAdded = 0;
	mat Q  = zeros(m, m );
	mat Qrho  = zeros(m, m );
	vec psi = zeros(m) ; // corresponding to each state over each interval
	vec psi_time = zeros(n); // for each tip, varies over time 
	
	vec A_Y;
	vec Y ;
	vec A = zeros(m);
	mat F = zeros(m,m);
	mat G = zeros(m,m);
	vec puY, pvY, pa;
	vec m_rs_R; 
	std::vector<bool> extant(n + Nnode, false); 
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
	DQAL dqal(Fs, Gs, Ys, m, hres, treeT ); 
	ode_state_type x0;//, x; 
	
	DQpsiA dqpsia( Fs, Gs, Ys, m, hres, treeT );
	ode_state_type x0sa; 
	
	// iterate events 
	double nextEventHeight = eventHeights(0); 
	int ievent = 0; 
	
	double sterm0, sterm1; 
	double h, h0, h1, hstar, dh; 
	int ih; // index of nextEventHeight
	double default_step_size = std::abs(heights[1] - heights[0]); 
	double L = 0.; 	
	
	//
	//~ boost::numeric::odeint::euler< ode_state_type > eu_stepper;
	
	h = 0.; 
	while( nextEventHeight != INFINITY && nextEventHeight < maxHeight ){
		if (nextEventHeight > h )
		{
			A = sum(P, 1);
			A = arma::normalise(A,1.)  * ((double)nextant);
			x0 = generate_initial_conditions(A) ; 
			size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight 
			  , std::min( default_step_size,  (nextEventHeight-h)/10.) );  
			
			//  same for sa ... 
			x0sa = generate_initial_conditions_dqpsia( A) ; 
			steps = boost::numeric::odeint::integrate( dqpsia ,  x0sa , h , nextEventHeight 
			  , std::min( default_step_size,  (nextEventHeight-h)/10.) );  
			Q_from_state(Q, x0); 
			A_from_state(A, x0); 
			L = L_from_state(x0); 
			A = arma::normalise(A,1.) * ((double)nextant);
			
			Qrho_from_state( Qrho, x0sa ); 
			psi_from_state( psi, x0sa );
			
			//update psi_time 
			// note need to run this before rho updated
			//~ double psi_factor; 
			//~ for (u = 0; u < n ; u++){
				//~ psi_factor = 0.; 
				//~ for (k = 0; k < m; k++){
					//~ psi_factor += rho.at(k,u) * psi.at(k); 
				//~ }
				//~ //psi_factor /= (double)m; 
				//~ psi_time.at(u) *= psi_factor;
			//~ }
			//vec psi_factor =  rho.t() * psi  ;
			psi_time = psi_time % (rho.t() * psi ); 
			
			P  = abs(Q.t() * P);
			P = arma::normalise( P, 1., 0 );
			rho = abs( Qrho.t() * rho ); 
			rho = arma::normalise( rho, 1., 0); // p=1, dim=0
		} else{
			L = 0.; 
		}
				
		if (eventIndicator(ievent)==SAMPLE)
		{
			u = eventIndicatorNode(ievent); 
			nextant++; 
			//~ P.col(u-1) = arma::normalise(sortedSampleStates.col(samplesAdded) ,1.); //TODO
			P.col(u-1) = arma::normalise(sampleStates.col(u-1) ,1.);
			rho.col(u-1) = P.col(u-1); 
			psi_time.at(u-1) = 1.; 
			extant.at(u-1) = true; 
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
			puY = arma::normalise(  arma::min(Y,P.col(u- 1 )) ,1.) / arma::clamp( Y, 1e-6, INFINITY ) ; 
			pvY = arma::normalise(  arma::min(Y,P.col(v- 1 )) ,1.) / arma::clamp( Y, 1e-6, INFINITY ) ; 
			//~ loglik += log( (( puY * F) * pvY).at(0,0) ) ; 
			
			// state of ancestor 
			pa =  arma::normalise( (F * puY) % pvY + (F * pvY) % puY ,1.) ; 
			//~ pa = pa / sum(pa ) ; 
			P.col(a - 1 ) = pa; 
			//P = finite_size_correction2(pa, A, extant, P);
			
			//bookkeeping
			nextant--; 
			extant.at(u-1) = false;
			extant.at(v-1) = false;
			extant.at(a-1) = true;
			P.col(u-1) = zeros<colvec>(P.n_rows); 
			P.col(v-1) = zeros<colvec>(P.n_rows);
			
			// sa stuff ; upate psi , rho 
			vec rho_w__Y = zeros(m); 
			vec rho_z__Y = zeros(m); 
			double pwz0,pzw0,pwz,pzw;
			// update W(iw, iz)
			for (int iw = 0; iw < n; iw++){
				if (tipsDescendedFrom.at(u-1, iw)==1){
					rho_w__Y = arma::normalise(  arma::min(Y,rho.col(iw)) ,1.) / arma::clamp(Y, 1e-16, INFINITY ) ; 
					for (int iz = iw+1; iz < n; iz++){
						if (tipsDescendedFrom.at(v-1, iz)==1){
							rho_z__Y = arma::normalise(  arma::min(Y,rho.col(iz)) ,1.) / arma::clamp(Y, 1e-16, INFINITY ) ; 
							pwz0 = sum( rho_w__Y % (F * rho_z__Y) );
							pzw0 = sum( rho_z__Y % (F * rho_w__Y ));
							
							if ( (pwz0 + pzw0) > 0 ){
								pwz = psi_time.at(iw) * psi_time.at(iz) * pwz0 / (pwz0 + pzw0) ;
								pzw = psi_time.at(iw) * psi_time.at(iz) *  pzw0 / (pwz0 + pzw0 ) ;
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
					psi_time.at(iw) *= puv ;
					//update rho(iw)
					rho.col(iw) = arma::normalise( rho_w__Y % (F * pvY) , 1.);
				} else if (tipsDescendedFrom.at(v-1, iw)==1){
					//tip descended from a and v
					rho_w__Y = arma::normalise(arma::clamp(  arma::min(Y,rho.col(iw))  / Y, 0., 1. ) ,1.) ; 
					//update psi(iw)
					double pvu = sum( rho_w__Y % (F * puY) );
					double puv = sum( puY % (F * rho_w__Y ));
					pvu = pvu / (puv + pvu ) ;
					psi_time.at(iw) *= pvu ;
					//update rho(iw)
					rho.col(iw) = arma::normalise(  rho_w__Y % (F * puY) ,1.);
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

