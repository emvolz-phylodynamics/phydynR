/* Helper functions for modular colik
 */

// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>

static const double MINY = 1e-12;

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
		if (u!=a){
			u = extantLines(iu); 
			p_u = mstates.col(u-1); 
			rterm = p_a / clamp(( A - p_u), 1e-12, INFINITY );
			rho = A / clamp(( A - p_u), 1e-12, INFINITY );
			lterm = dot( rho, p_a); //
			p_u = p_u % clamp((lterm - rterm), 0., INFINITY) ;
			p_u = p_u / sum(p_u ) ; 
			mstates.col(u - 1) = p_u; 
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
	vec pu__Y = clamp( pu/Y , 0., 1.);
	vec pv__Y = clamp( pv/Y , 0., 1.);
//~ cout << pu << endl;
//~ cout << pv << endl;
//~ cout << Y << endl;
//~ cout << pu__Y << endl;
//~ cout << pv__Y << endl;
	vec pa = pu__Y % (F * pv__Y) + pv__Y % (F*pu__Y); 
	//~ vec pa = pu__Y % dot(F , pv__Y) + pv__Y % dot(F, pu__Y); //dnw
	double corate = sum(pa); 
	pa = normalise( pa, 1.);
//~ cout << pa << endl;
//~ cout << "$$$$$$$$$$$$" << endl;
	List o; 
	o["pa"] = pa;
	o["corate"] = corate;
	return o; 
}


//~ .update.alpha <- function(u,v, tree, fgy, A)
//~ { 
	//~ .F <- fgy$.F
	//~ .G <- fgy$.G
	//~ .Y <- fgy$.Y
	//~ {
		//~ .Y <- pmax(A, .Y)
		//~ FklXpuk_Yk <- (.F * tree$mstates[,u]/.Y)
		//~ FklXpvk_Yk <- (.F * tree$mstates[,v]/.Y)
		//~ FklXpuk_Yk[is.nan(FklXpuk_Yk)] <- 0
		//~ FklXpvk_Yk[is.nan(FklXpvk_Yk)] <- 0
		//~ vk_Yk <- pmin(pmax(tree$mstates[,v]/.Y, 0),1); vk_Yk[is.nan(vk_Yk)] <- 0
		//~ uk_Yk <- pmin(pmax(tree$mstates[,u]/.Y, 0),1); uk_Yk[is.nan(uk_Yk)] <- 0
		//~ ratekl <- FklXpuk_Yk %*% vk_Yk + FklXpvk_Yk %*% uk_Yk
	//~ }
	//~ 
	//~ coalescentRate <- max( sum(ratekl) , 0)
	//~ 
	//~ if (sum(ratekl)==0) {ratekl <- rep(1/tree$m, tree$m) * 1e-6}
	//~ # definitions of alpha state
	//~ p_a <- as.vector( ratekl / sum(ratekl) )
//~ #~ if ( tree$heights[u] > 11) browser()
	//~ list( pa = p_a, corate =coalescentRate) 
//~ }
//~ 


//[[Rcpp::export()]]
void update_states0( mat& mstates , mat Q)
{
	// mstates m x (n+Nnode)
	mstates = Q.t() * mstates; 
	mstates = normalise( clamp(mstates, 0., 1.), 1., 0);
}

