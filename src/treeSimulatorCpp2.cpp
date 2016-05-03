/* 
 * rcolgem tree simulator 
 * 
 */


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::plugins(cpp11)]]



static const double A_EPS = 0.00;
static const double Y_EPS = 1.E-12;


using namespace std;
using namespace Rcpp; 
using namespace arma; 


mat finite_size_correction2(const vec& p_a, const vec& A, const IntegerVector extant, mat mstates)
{ //compare to version in colik.cpp
// TODO not working yet
	// NOTE mstates m X n 
	vec rho; 
	vec rterm; 
	//~ vec lterm; 
	double lterm; 
	vec p_u; 
	vec Amin_pu; 
	//~ for (int iu = 0; iu < extantLines.size(); iu++){
	for (int iu = 0; iu < extant.size(); iu++){
		if (extant(iu)==1){
			p_u = mstates.row(iu); 
			Amin_pu = clamp(( A - p_u), 1., INFINITY ); 
			rterm = p_a / Amin_pu ; 
			rho = A / Amin_pu; 
			lterm = dot( rho, p_a); //
			p_u = p_u % (lterm - rterm) ;
			p_u = p_u / sum(p_u ) ; 
			mstates.row(iu) = p_u; 
		}
	}
	return mstates; 
}


void incorporateSamples( double h, int samplesAdded, mat lstates, mat mstates, NumericVector sortedSampleHeights, NumericMatrix sortedSampleStates, IntegerVector extant,  NumericVector heights, int m ){
	for (int i = samplesAdded; i < sortedSampleHeights.size(); i++){
		if (sortedSampleHeights[i]==h){
			extant[i] = 1; 
			heights[i]  =h;
			for (int k = 0; k < m; k++) {
				lstates(i,k) = sortedSampleStates(i,k);
				mstates(i,k) = sortedSampleStates(i,k);
			}
		}
		if (sortedSampleHeights[i] > h) break;
	}
}

//[[Rcpp::export]]
List simulateTreeCpp2(const NumericVector times, const List Fs, const List Gs, const List Ys
  , const List As
  , NumericVector sortedCoHeights // computed at R level
  , const NumericVector sortedSampleHeights
  , const NumericMatrix sortedSampleStates
  , double maxSampleTime 
  , const int m
  , bool finiteSizeCorrection
  , std::vector< std::string > DEMES
)
{
	RNGScope scope;
	
	int n = sortedSampleHeights.size();
	int Nnode = sortedCoHeights.size(); //n-1; 
	// note times decreasing
	double tfin = times[times.size()-1];
	double maxHeight = maxSampleTime -tfin;
	const double deltat  = times[0] - times[1]; 
	
	IntegerVector extant( n + Nnode, 0);
	IntegerVector extantRange(n+Nnode);
	for  (int i = 0; i < n+Nnode; i++) extantRange[i] = i; 
	int numberExtant = 0;
	
	// note root does not have edge
	int ntrees = 1 + ((n -1) - Nnode); 
	NumericMatrix edge(n+Nnode-ntrees, 2); 
	NumericVector edge_length(n+Nnode-ntrees, -1.0);
	NumericVector heights(n+Nnode, -1.0);
	
	mat lstates = zeros(n+Nnode, m);
	mat mstates = zeros(n+Nnode, m);
	mat ustates = -1. * ones(n+Nnode, m);
	IntegerVector donor( n + Nnode, -1 ); // which daughter node is donor at each internal node
	IntegerVector recip( n + Nnode , -1); 
	
	double t, t0, t1;
	int internalNodeIndex = n; // counter for internal branches
	int internalNodesAdded = 0;
	int samplesAdded = 0; // counter for terminal branches
	int edgesAdded = 0;
	
	vec A_Y;
	vec Y ;
	vec A;
	mat F = zeros(m,m);
	mat G = zeros(m,m);
	int it, ist, k, l, nco, inco , i, stateOf_a, stateOf_recip;
	int u, v, a, z, w;
	int nextant; 
	double r, rk, rl, cwk, cwl, twk, twl, h1;
	bool foundv, foundu;
	double rsmsi;
	
	double h;
		
	h = 0.; 
	
	for (int i = samplesAdded; i < sortedSampleHeights.size(); i++){
		if (sortedSampleHeights[i]==h){
			extant[i] = 1; 
			heights[i]  =h;
			for (int k = 0; k < m; k++) {
				lstates(i,k) = sortedSampleStates(i,k);
				mstates(i,k) = sortedSampleStates(i,k);
			}
			samplesAdded++;
			numberExtant++; 
		}
		if (sortedSampleHeights[i] > h) break;
	}
	
	const double deltah  = times[0] - times[1]; //note times decreasing
	
	NumericMatrix corates(m,m);
	mat R  = zeros(m, m );
	mat Q  = zeros(m, m );
	vec m_rs_R;
	
	// note times decreasing
	for(it = 0; it < times.size()-1; it++)
	{
		t =  times[it];
		t1 =times[it+1];  //t - deltat;
		h = maxSampleTime - t; 
		h1 = maxSampleTime - t1; 

		
		F = as<mat>(Fs[it]);
		G = as<mat>(Gs[it]);
		Y = as<vec>(Ys[it]);
		A = as<vec>(As[it]);
		// transition probabilities 
		R = trans( F + G ) ;
		for ( k = 0; k < Y.size(); k++)
		{
			//~ if (Y.at(k) > 0){
				R.col(k) = R.col(k) / clamp(Y, 1e-6, std::max(1e-6, Y.max()));
			//~ } else{
				//~ R.col(k).fill(0.);
				//~ R.col(k) = zeros<colvec>(Y.size());
			//~ }
		}
		R.diag(0) =   zeros<colvec>(Y.size());
		m_rs_R = -sum( R, 1 );
		// diagonal correction for coalescent events
		A_Y = clamp( A / clamp(Y, 1e-6, std::max(1e-6,Y.max())), 0., 1. ); 
		vec diag_xn = zeros<colvec>(m); 
		//~ for (k = 0; k < m; k++){
			//~ diag_xn[k] = ((double)dot(F.row(k) , A_Y))/ std::max(Y.at(k), 1.) ;//Y[k]; 
		//~ }
		R.diag(0) = m_rs_R; //- diag_xn; //TODO 
		// transition prob from row to col
		Q = expmat( (h1 - h) * R); 
		// renormalise
		for (k = 0; k < m; k++){
			double rsq = sum(Q.row(k));
			if (rsq > 0){
				Q.row(k) = Q.row(k) / rsq;
			}else{
				Q(k,k) = 1. ;
			}
		}
		// </transition probs>
		
		//update line states
		for ( i = 0; i < (n + internalNodesAdded); i++)
		{
			if (0!=extant[i]){
				mstates.row(i) = mstates.row(i) * Q; 
				//~ for ( k = 0; k < m; k++)
				//~ {
					//~ mstates(i,k) = std::max(0.0, (double)(dot(Q.col(k), mstates.row(i))) );
				//~ }
				//renormalise
				rsmsi  = sum(mstates.row(i));
				for ( k =0; k < m; k++){
					mstates(i,k) /= rsmsi;
				}
			}
		}
		
		
		// co events
		double corate = 0.;
		for ( k = 0; k < m; k++) {
			//~ corates.row(k) =  A_Y[k] * (F.row(k) % A_Y)
			for ( l = 0; l < m; l++){
				if (k == l){
					corates(k,k) =  A_Y(k) * (std::max(1.0,(A[k]-1))/std::max(1.0,(Y[k]-1)))  * F(k,k);
				} else{
					corates(k,l) =  A_Y(k) * A_Y(l) * F(k,l);
				}
				corate += corates(k,l);
			}
		}
		//~ for (i = internalNodesAdded; i < sortedCoHeights.size(); i++){
		for (int ico = 0; ico < sortedCoHeights.size(); ico++){
			if (sortedCoHeights[ico] > h1) break;
			// ensure that at least two lines are extant if doing co: 
			nextant = sum(extant); 
			if ( nextant < 2 && (sortedCoHeights[ico] > h && sortedCoHeights[ico] <= h1) ){
				sortedCoHeights[ico] = h1+(h1-h)*1e-3; 
			}
			if (sortedCoHeights[ico] > h && sortedCoHeights[ico] <= h1){
				// find k and l
				r = R::runif( 0., corate); 
				double sumcorate = 0.;
				int donordeme = 0;
				int recipdeme = 0;
				bool demes_found = false;
				for (w = 0; w < m; w++){
					for (z = 0; z < m; z++){
						sumcorate += corates(w,z);
						if (sumcorate > r){
							donordeme = w;
							recipdeme = z;
							demes_found = true; 
							break;
						}
					}
					if (sumcorate > r){
						break;
					}
				}
				// sample u & v in proportion to mstates k & l
				foundu = false;
				foundv = false;
				double Ak = 0.;//sum(mstates(_,k));
				double Al = 0.;//sum(mstates(_,l));
				for (i = 0; i < (n+internalNodesAdded); i++){
					if (extant[i]!=0){
						Ak += mstates(i, donordeme);
						Al += mstates(i,recipdeme);
					}
				}
				rk = R::runif(0., Ak); 
				double sumAk =0.;
				double sumAl = 0.;
				for (i = 0; i < (n + internalNodesAdded); i++){
					if (extant[i]!=0){
						sumAk += mstates(i,donordeme);
						if (sumAk > rk){
							u  = i; 
							foundu = true;
							Al -= mstates(i,recipdeme);
							break;
						} 
					}
				}
				rl = R::runif( 0., Al); 
				for (i = 0; i < (n + internalNodesAdded); i++){
					if ((i != u) && (extant[i]!=0)){
						sumAl += mstates(i,recipdeme);
						if (sumAl > rl){
							v = i;
							foundv = true;
							break; 
						} 
					}
				}
				
				// do the co
				if (!foundu || !foundv){
					// just pick two lines at random; 
					//~ cout << sum(extant) << endl; 
					std::cout << "Warning: could not find compatible pair of lines for coalescent" << std::endl; 
					IntegerVector uv =  Rcpp::RcppArmadillo::sample( extantRange , 2, false, as<NumericVector>(extant)) ;
					u = uv[0];
					v = uv[1];
				}
//~ std::cout << donordeme << " " << recipdeme << std::endl; 
//~ std::cout << DEMES.at(donordeme) << " "<< DEMES.at(recipdeme) << std::endl; 
//~ std::cout << internalNodesAdded << " " << samplesAdded << std::endl; 
				a = n + internalNodesAdded;
				internalNodesAdded++;
				heights[a] = sortedCoHeights[ico]; 
				
				edge(edgesAdded, 0) = a;
				edge(edgesAdded, 1) = u;
				edge_length[edgesAdded] = heights[a]-heights[u];
				edgesAdded++;
				edge(edgesAdded, 0) = a;
				edge(edgesAdded, 1) = v;
				edge_length[edgesAdded] = heights[a]-heights[v];
				edgesAdded++;
				
				extant[u] = 0;
				extant[v] = 0;
				extant[a] = 1;
				//~ numberExtant--;
				
				//lstate, mstate of a;
				/* TODO does  this gives overly ladder-like trees??
				lstates(a,donordeme) = 1.0;
				mstates(a,donordeme) = 1.0;
				*/
				for ( int w = 0; w < m ; w++){
					for (int z = 0; z < m; z++){
						lstates(a, w) = F(w,z) * ( mstates(u,w)*mstates(v,z) / Y(z)/Y(w) +  mstates(u,z)*mstates(v,w) / Y(z)/Y(w) ) ;
					}
				}
				lstates.row(a) = lstates.row(a) / sum(lstates.row(a)); 
				mstates.row(a) = lstates.row(a); 
				
				// set ustates for u and v
				ustates.row(u) = mstates.row(u); 
				ustates.row(v) = mstates.row(v); 
				donor(a) = u; 
				recip(a) = v; 
				
				//update mstates of lines not in co 
				if (finiteSizeCorrection){
					std::cout << "finiteSizeCorrection not implemented" << std::endl; 
					throw 1; 
					mstates = finite_size_correction2(lstates.row(a), A, extant, mstates);
					/*
					for (int s = 0; s < internalNodesAdded-1; s++){
						if (extant[s]!=0){
							for (int w = 0; w <m; w++){
								if (w!=k ){
									if (mstates(s,k) > 0){
										mstates(s, w) = mstates(s,w) * (Y[k] / (Y[k] - mstates(s,k)) )  ;
									}
								}
								if (w!=l ){
									if (mstates(s,l) > 0){
										mstates(s, w) = mstates(s,w) * (Y[l] / (Y[l] - mstates(s,l)) )  ;
									}
								}
							}
							mstates(s,k) = mstates(s,k) * std::max(0., ((Y[k]-1.)/Y[k])) ;
							mstates(s,l) = mstates(s,l) * std::max(0., ((Y[l]-1.)/Y[l])) ;
							mstates.row(s) = mstates.row(s) / sum(mstates.row(s));
						}
					}
					*/
				}
			}
		}
		
		// terminate early if all nodes added 
		if ( internalNodesAdded >= sortedCoHeights.size() ) break; 
		
		// find sample times in (t, t1)
		//~ incorporateSamples(   h,  samplesAdded,  lstates,  mstates,  sortedSampleHeights,  sortedSampleStates,  extant,   heights, m);
		// new samples
		for (int i = samplesAdded; i < sortedSampleHeights.size(); i++){
			if ((sortedSampleHeights[i]>h) && (sortedSampleHeights[i]<=h1)){
				extant[i] = 1; 
				heights[i]  = sortedSampleHeights[i];
				for (int k = 0; k < m; k++) {
					lstates(i,k) = sortedSampleStates(i,k);
					mstates(i,k) = sortedSampleStates(i,k);
				}
				samplesAdded++;
				numberExtant++;
			}
			if (sortedSampleHeights[i] > h1) break;
		}
		
		
	}
	
	List ret;
	ret["edge"] = edge;
	ret["edge.length"] = edge_length;
	ret["n"] = n;
	ret["Nnode"] = internalNodesAdded;//Nnode
	ret["lstates"] = lstates;
	ret["mstates"] = mstates;
	ret["ustates"] = ustates;
	ret["donor"] = donor;
	ret["recip"] = recip;
	ret["heights"] = heights;
	ret["samplesAdded"] = samplesAdded; 
	
	// debug output
	ret["R"]  = R; 
	ret["Q"] = Q; 
	//~ ret["mstates"] = mstates; 
	
	return(ret);
}
