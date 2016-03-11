/* lineages through time and hazard of co
 */

#include "colik.hpp" 

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, BH)]]


using namespace std;
using namespace Rcpp; 
using namespace arma;

//~ static const int SAMPLE  = 0; 
//~ static const int CO = 1; 



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


// transition probabilities and survivor functions, conditioning on lineage not changing host
class DQpsiA{
	List Fs, Gs, Ys; 
	int m; 
	double hres; 
	double treeT; 
public:
	DQAL( List Fs, List  Gs, List Ys, int m, double hres, double treeT ) : Fs(Fs),Gs(Gs),Ys(Ys),m(m),hres(hres),treeT(treeT) {};
	void operator() ( const state_type &x , state_type &dxdt , double t)//, const double /* t */ )
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
				for (l = 0. ; l < m; l++){
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
		// d psi  
		for (z = 0; z < m; z++){ // col of psi- line starts in deme z
			for (k = 0; k < m; k++){ //row of psi - current deme
				dxdt[ psi_ind(k, z ) ] = 0. ; 
				for (l = 0. ; l < m; l++){ //source of transm
					// TODO possibly should have (Y_l - A_l)/Y_l terms as in paper 
					//~ dxdt[psi_ind(k,z)] -= F(i,l,k) * psi(x, k, z) / std::max(psi(x, k, z), Y(k)); 
					// NOTE Q is filling the role of rho_{ik} in the paper 
					dxdt[psi_ind(z)] -= F(i,l,k) * psi(x, z) * Q(x,k,z) / std::max(Q(x, k, z), Y(k)); 
				}
			}
		}
		// NO
		//~ for (k = 0; k < m; k++){
			//~ for (l = 0; l < m; l++){
				//~ dxdt[psi_ind( l )] -= F(i,l,k) * psi(x, k) / std::max(psi(x,k), Y(k)); 
			//~ }
		//~ }
    }    
private:
	double Q( const state_type &x, int k, int l) {
		return x[l*m+k]; 
	}
	double psi( const state_type &x, int k) {
		return x[(int)pow(m,2) +  k]; 
	}
	double A( const state_type &x, int k) {
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



state_type generate_initial_conditions_dqpsia(vec A){
	int m = A.size() ;
	state_type x( (int)pow(m,2) + m + m, 0. ); 
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


void Qrho_from_state( mat &Qrho, state_type xfin ){
	int m = Q.n_rows; 
	int k =0 ;
	for (int i = 0; i < m; i++){
		for (int j = 0; j < m; j++){
			Qrho.at(i,j) = xfin[k]; 
			k++; 
		}
	}
	for (int i = 0; i < m; i++){
		Q.row(i) = Q.row(i) / sum(Q.row(i)); 
	}
}

void psi_from_state( vec &psi, state_type xfin ){
	int m = psi.size(); 
	int k = 0; 
	for ( int i = (int)pow(m,2); i < ((int)pow(m,2)+m); i++){
		psi.at(k) = xfin[i]; 
		k++; 
	}
}

//[[Rcpp::export()]]
mat sourceAttribMultiDemeCpp(const NumericVector heights, const List Fs, const List Gs, const List Ys
  , const IntegerVector eventIndicator // sample or co
  , const IntegerVector eventIndicatorNode // node involved at each event
  , const NumericVector eventHeights
  , const mat sortedSampleStates
  , const IntegerMatrix daughters // daughters of each node
  , const int n
  , const int Nnode
  , const int m
  , bool AgtYboundaryCondition
  , const double maxHeight // terminate computation at this height
)
{
	double loglik = 0.;
	mat P(m, n + Nnode, fill::zeros); 
	mat rho(m, n , fill::zeros); 
	int nextant = 0; 
	int u,v,w,z,a, k, l; 
	int samplesAdded = 0;
	mat Q  = zeros(m, m );
	mat Qrho  = zeros(m, m );
	//~ vec Psi = zeros(n, n + Nnode); //records probability i'th sample is host at j'th node 
	vec psi = zeros(m) ; // corresponding to each state over each interval
	vec psi_time = zeros(n); // for each tip, varies over time 
	mat W = zeros(n,n); 
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
	
	// instantiate solver 
	double hres =  heights.size() ;
	double treeT = heights[heights.size()-1]; // note increasing order   
	DQAL dqal(Fs, Gs, Ys, m, hres, treeT ); 
	state_type x0;//, x; 
	
	DQpsiA dqpsia( Fs, Gs, Ys, m, hres, treeT );
	state_type x0sa; 
	
	// iterate events 
	double nextEventHeight = eventHeights(0); 
	int ievent = 0; 
	
	double sterm0, sterm1; 
	double h, h0, h1, hstar, dh; 
	int ih; // index of nextEventHeight
	double default_step_size = std::abs(heights[1] - heights[0]); 
	double L = 0.; 	
	
	h = 0.; 
	while( nextEventHeight != INFINITY ){
		
		if (nextEventHeight > h ){
			A = sum(P, 1);
			A = normalise(A)  * ((double)nextant);
			x0 = generate_initial_conditions(A) ; 
//~ cout << "solving..." << endl;
//~ cout << h << endl;
//~ cout << nextEventHeight << endl;
//~ cout << A << endl;
//~ Rf_PrintValue( wrap(x0)); 
			//~ size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight , default_step_size, simple_observer(x)  );  
			size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight 
			  , std::min( default_step_size,  (nextEventHeight-h)/10.) );  
			//~ size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight , default_step_size
			  //~ , push_back_state_and_time( x_vec, x_vec_times)  );  
			//~ size_t steps = boost::numeric::odeint::integrate( dqal ,  x0 , h , nextEventHeight , (nextEventHeight-h)/10.
			  //~ , push_back_state_and_time( x_vec, x_vec_times)  );  
//~ Rf_PrintValue( wrap(x0)); 
//~ cout << steps << endl;
//~ throw 1; 
			// 
			
			//  same for sa ... 
			x0sa = generate_initial_conditions_dqpsia( A) ; 
			size_t steps = boost::numeric::odeint::integrate( dqpsia ,  x0sa , h , nextEventHeight 
			  , std::min( default_step_size,  (nextEventHeight-h)/10.) );  
			
			Q_from_state(Q, x0); 
			A_from_state(A, x0); 
			L = L_from_state(x0); 
			A = normalise(A) * ((double)nextant);
			
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
			rho = abs( Qrho.t() * rho ); 
			
			
		} else{
			L = 0.; 
		}
				
		if (eventIndicator(ievent)==SAMPLE)
		{
			u = eventIndicatorNode(ievent); 
			nextant++; 
			P.col(u-1) = normalise(sortedSampleStates.col(samplesAdded));
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
//~ if (sum(A) > sum(Y)) {
	//~ cout << Y ; 
	//~ cout << nextEventHeight << " " << ih << endl; 
	//~ cout << " sum A > sum Y  " << endl ;
//~ }
			
			a = eventIndicatorNode(ievent);
			u = daughters(a -1 , 0); 
			v = daughters( a - 1, 1 ); 
			puY = normalise(  arma::min(Y,P.col(u- 1 )) ) / arma::clamp(Y, 1e-6, INFINITY ) ; 
			pvY = normalise(  arma::min(Y,P.col(v- 1 )) ) / arma::clamp( Y, 1e-6, INFINITY ) ; 
			//~ loglik += log( (( puY * F) * pvY).at(0,0) ) ; 
			
			// state of ancestor 
			pa =  arma::normalise( (F * puY) % pvY + (F * pvY) % puY) ; 
			//~ pa = pa / sum(pa ) ; 
			P.col(a - 1 ) = pa; 
			P = finite_size_correction2(pa, A, extant, P);
			
			//bookkeeping
			nextant--; 
			extant.at(u-1) = false;
			extant.at(v-1) = false;
			extant.at(a-1) = true;
			P.col(u-1) = zeros<colvec>(P.n_rows); 
			P.col(v-1) = zeros<colvec>(P.n_rows);
			
			// sa stuff ; upate psi , rho 
			tipsDescendedFrom.row(a-1) = tipsDescendedFrom.row(u-1) + tipsDescendedFrom.row(v-1); 
			for (int iw = 0; iw < n; iw++){// w & z tips
				if (tipsDescendedFrom.at(a-1, iw)==1){
					w= iw + 1 ; 
					//update psi(iw)
					
					//update rho(iw)
				}
			}
			for (int iw = 0; iw < n; iw++){
				if (tipsDescendedFrom.at(a-1, iw)==1){
					for (int iz = iw+1; iz < n; iz++){
						if (tipsDescendedFrom.at(a-1, iz)==1){
							// update W(iw, iz)
						}
					}
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
	}
	
	return W; 
}


