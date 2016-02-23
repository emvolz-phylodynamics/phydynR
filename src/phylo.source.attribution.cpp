/* lineages through time and hazard of co
 */
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]


using namespace std;
using namespace Rcpp; 

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
