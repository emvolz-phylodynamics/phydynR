/* Helper functions for modular colik
 */

// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>

//~ static const double MINY = 1e-12;
static const double MINY = 1.;

using namespace arma;
using namespace Rcpp; 
using namespace std; 

//[[Rcpp::export]]
void rco_finite_size_correction2(const int a, const vec& p_a, const vec& A, const uvec& extantLines, mat& mstates)
{
	// NOTE mstates m X n 
	int u; 
	vec rho; 
	vec rterm; 
	//~ vec lterm; 
	double lterm; 
	vec p_u; 
	for (int iu = 0; iu < extantLines.size(); iu++){
		u = extantLines(iu); 
		if (u!=a){
			p_u = mstates.col(u-1); 
			rterm = p_a / clamp(( A - p_u), 1e-12, INFINITY );
			rho = A / clamp(( A - p_u), 1e-12, INFINITY );
			lterm = dot( rho, p_a); //
			p_u = p_u % clamp((lterm - rterm), 0., INFINITY) ;
			if (sum(p_u) > 0.){
				p_u = p_u / sum(p_u ) ; 
				mstates.col(u - 1) = p_u; 
			}
		}
	}
}


//[[Rcpp::export()]]
List eventTimes2extant( vec eventTimes, vec nodeheights, vec parentheights )
{
	List o(eventTimes.size());
	List nodesAtHeight( eventTimes.size());  
	double h, nodeheight, parentheight ; 
	for (int i =0  ; i < eventTimes.size(); i++){
		h = eventTimes(i); 
		std::vector<int> x; 
		x.reserve(nodeheights.size()); 
		std::vector<int> nah; 
		nah.reserve( nodeheights.size()); 
		for (int u = 0 ; u < nodeheights.size(); u++){
			nodeheight = nodeheights(u);
			parentheight = parentheights(u); 
			if (nodeheight==h){
				nah.push_back( u + 1 );
			}
			if (!NumericVector::is_na( parentheight)){
				if ( nodeheight<= h && parentheight > h ){
					x.push_back( u+1 ); 
				}	
			} else{
				if (nodeheight <= h ){
					x.push_back( u+1 ) ; 
				}
			}
		}
		o[i] = wrap(x); 
		nodesAtHeight[i] = wrap(nah); 
	} 
	List oo; 
	oo["extant"] = o;
	oo["nodesAtHeight"] = nodesAtHeight ; 
	return oo;
}


// NOTE also tallies rate and survival likelihood terms
//[[Rcpp::export()]]
List update_alpha0(vec pu
  , vec pv
  , mat F
  , vec Y
  , vec A
){
	int m = Y.size(); 
	Y = max( A, clamp(Y, MINY, INFINITY)); //NOTE this line is very important! 
	//~ vec pu__Y = clamp( pu/Y , 0., 1.);
	//~ vec pv__Y = clamp( pv/Y , 0., 1.);
	vec pu__Y =pu/Y ;
	vec pv__Y =pv/Y ;
	vec pa = pu__Y % (F * pv__Y) + pv__Y % (F*pu__Y); 
	double corate = sum(pa); 
	//~ double corate =  F(0,0)   / Y(0) / Y(0) / 2.; 
	pa = normalise( pa, 1.);
	List o; 
	o["pa"] = pa;
	o["corate"] = corate;
	return o; 
}


//[[Rcpp::export()]]
void update_states0( mat& mstates , mat Q)
{
	// mstates m x (n+Nnode)
	mstates = Q.t() * mstates; 
	mstates = normalise( clamp(mstates, 0., 1.), 1., 0);
}

//[[Rcpp::export()]]
mat update_states1( mat& mstates, mat Q , vec extantLines){
	int u; 
	for (int iu = 0; iu < extantLines.size(); iu++){
		u = extantLines(iu)-1; // R -> c indexing
		mstates.col(u)  = normalise( clamp( Q.t() * mstates.col(u), 0., 1.), 1.);
	}
	return mstates;
}
