/* Helper functions for modular colik
 */

// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>

static const double MINY = 1e-12;

using namespace arma;
using namespace Rcpp; 


//[[Rcpp::export()]]
void finite_size_correction0(mat& mstates
  , int a
  , vec pa
  , uvec extantLines
  , vec A
){ // mods mstates
	int ia = a - 1;
	//update P
	// NOTE mstates m X n 
	vec rho; 
	vec rterm; 
	//~ vec lterm; 
	double lterm; 
	vec p_u; 
	vec Amin_pu; 
	for (int iiu = 0; iiu < extantLines.size(); iiu++){
		int iu = extantLines(iiu); 
		p_u = mstates.row(iu); // TODO to colvec? 
		//~ Amin_pu = clamp(( A - p_u), 1., INFINITY ); 
		Amin_pu = clamp(( A - p_u), 1e-3, INFINITY ); 
		rterm = pa / Amin_pu ; 
		rho = A / Amin_pu; 
		lterm = dot( rho, pa); //
		p_u = p_u % clamp((lterm - rterm), 0., INFINITY) ; // l > r
		p_u = p_u / sum(p_u ) ; 
		//~ P.col(iu) = p_u; 
		mstates.row(iu) = p_u; // TODO to row?
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
				if ( nodeheight>= h && parentheight > h ){
					x.push_back( u+1 ); 
				}	
			} else{
				if (nodeheight >= h ){
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
){
	int m = Y.size(); 
	Y = clamp(Y, MINY, INFINITY); 
	vec pu__Y = clamp( pu/Y , 0., 1.);
	vec pv__Y = clamp( pv/Y , 0., 1.);
	vec pa = pu__Y % (F * pv__Y) + pv__Y % (F*pu__Y); 
	double corate = sum(pa); 
	pa = normalise( pa, 1.);
	List o; 
	o["pa"] = pa;
	o["corate"] = corate;
	return o; 
}

//[[Rcpp::export()]]
void update_states0( mat& mstates , mat Q)
{
	mstates = Q.t() * mstates; 
	mstates = normalise( clamp(mstates, 0., 1.), 1., 0);
}

