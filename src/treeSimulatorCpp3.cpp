/* 
 * rcolgem tree simulator 
 * 
 */


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]



static const double A_EPS = 0.00;
static const double Y_EPS = 1.E-12;


using namespace std;
using namespace Rcpp; 
using namespace arma; 

class CoSim{
public://TODO after debugging make these private
	int n ;
	int Nnode ; 
	// note times decreasing
	double tfin;
	double maxHeight;
	double deltat; 
	
	vector<bool> extant;
	int i, j;
	int numberExtant;
	
	// note root does not have edge
	mat edge; 
	vec edge_length;
	vec heights;
	
	// for sim genetic distances
	vec edge_length_substPerSite; 
	vec edge_length_substPerSite_nodeWise;
	
	mat lambda_uv;
	
	mat P;
	mat lstates;
	mat ustates;
	
	double t, t0, t1;
	int internalNodeIndex ; // counter for internal branches
	int internalNodesAdded;
	int samplesAdded ; // counter for terminal branches
	int edgesAdded ;
	
	vec A;
	vec A_Y;
	vec Y ;
	mat F; 
	mat G; 
	int it, ist, k, l, nco, inco ;
	int u, v, a, z, w;
	double r;
	
	double h, h0, h1;
	
	//inputs
	vec times;
	List Fs, Gs,Ys; //TODO need to change? make global?
	vec sortedSampleHeights;
	mat sortedSampleStates; // m X n 
	double maxSampleTime ;
	int m;
	bool finiteSizeCorrection;
	vec substitutionRates; // m rates 
	int sequenceLength; 
	
	CoSim( vec in_times, List in_Fs,  List in_Gs,  List in_Ys
	  ,  vec in_sortedSampleHeights
	  ,  mat in_sortedSampleStates // m X n 
	  , double in_maxSampleTime 
	  , int in_m
	  , bool in_finiteSizeCorrection
	  , vec in_substitutionRates
	  , int in_sequenceLength) : times(in_times),Fs(in_Fs),Gs(in_Gs),Ys(in_Ys),sortedSampleHeights(in_sortedSampleHeights),sortedSampleStates(in_sortedSampleStates),maxSampleTime(in_maxSampleTime),m(in_m),finiteSizeCorrection(in_finiteSizeCorrection),substitutionRates(in_substitutionRates),sequenceLength(in_sequenceLength)
	{
		
		n = sortedSampleHeights.size();
		Nnode = n -1; 
		// note times decreasing
		tfin = times(times.size()-1);
		maxHeight = maxSampleTime -tfin;
		
		//~ extant(n + Nnode, false); 
		extant.resize(n + Nnode); std::fill( extant.begin(), extant.end(), false); 
		numberExtant = 0;
		
		// note root does not have edge
		edge = zeros(n+Nnode, 2);  //
		std::fill( edge.begin(), edge.end(), NumericVector::get_na() ) ;
		edge_length = vec(n+Nnode-1); //NOTE no root edge
		std::fill( edge_length.begin(), edge_length.end(), NumericVector::get_na() ) ;
		edge_length_substPerSite = vec(n+Nnode-1); //
		std::fill( edge_length_substPerSite.begin(), edge_length_substPerSite.end(), 0. ) ;
		edge_length_substPerSite_nodeWise = vec(n + Nnode);
		std::fill( edge_length_substPerSite_nodeWise.begin(), edge_length_substPerSite_nodeWise.end(), 0. ) ;
		//~ heights = vec(n+Nnode); heights.fill(INFINITY);
		heights = zeros(n+Nnode); heights.fill(INFINITY);
		
		lambda_uv = zeros(n + Nnode, n + Nnode ); 
		P = zeros( m, n + Nnode );
		lstates = zeros(m, n + Nnode);
		ustates = zeros(m, n + Nnode);
		
		internalNodeIndex = n; // counter for internal branches
		internalNodesAdded = 0;
		samplesAdded = 0; // counter for terminal branches
		edgesAdded = 0;
		
		A = zeros(m); 
		F = zeros(m,m);
		G = zeros(m,m);
		
		h = 0.; 
	}
	
	void incorporate_samples(const double h0
		 , const double h1
	)
	{
		for (int i = samplesAdded; i < sortedSampleHeights.size(); i++){
			if ( (sortedSampleHeights(i)>=h0) && (sortedSampleHeights(i)< h1) ){
				extant.at(i) = true; 
				heights(i)  =sortedSampleHeights(i);
				P.col(i) = sortedSampleStates.col(i); 
				lstates.col(i) = sortedSampleStates.col(i);
				samplesAdded++;
				numberExtant++;
			}
			if (sortedSampleHeights(i) > h1) {
				break;
			}
		}
	}
	
	
	void finite_size_correction2(int ia)
	{ //compare to version in colik.cpp
		//update P
		// NOTE mstates m X n 
		vec rho; 
		vec rterm; 
		//~ vec lterm; 
		double lterm; 
		vec p_u; 
		vec Amin_pu; 
		//~ for (int iu = 0; iu < extantLines.size(); iu++){
		for (int iu = 0; iu < extant.size(); iu++){
			if (extant.at(iu) && iu != ia){
				p_u = P.col(iu); 
				Amin_pu = clamp(( A - p_u), 1., INFINITY ); 
				rterm = P.col(ia) / Amin_pu ; 
				rho = A / Amin_pu; 
				lterm = dot( rho, P.col(ia)); //
				p_u = p_u % clamp((lterm - rterm), 0., INFINITY) ; // l > r
				p_u = p_u / sum(p_u ) ; 
				P.col(iu) = p_u; 
			}
		}
	}

	double co_rate_uv(const int u, const int v)
	{
		vec pu__Y = clamp(P.col(u),0., 1.)/clamp(Y, 1e-12, INFINITY) ;
		vec pv__Y = clamp(P.col(v), 0., 1.)/clamp(Y, 1e-12, INFINITY);
		//~ return as_scalar((pu__Y * (F * pv__Y) + pv__Y * (F * pu__Y))); 
		return as_scalar(dot(pu__Y , (F * pv__Y)) + dot(pv__Y , (F * pu__Y))); 
	}
	
	
	void update_states( const double dh){
		vec A_Y = clamp(A / Y, 0., 1.); 
		
		// make R & Q 
		//~ (F(k,l) + G(k,l)) *  Q(x, l,z)/  std::max(Q(x,l,z), Y(l));
		mat R = (F+G).t();
		//~ R  = R.each_col() / Y; 
		R.each_col() /= Y; 
		R.diag().zeros();
		R.diag() = -sum(R,1); //- ((F *  A_Y) / Y) ;
		mat Q = normalise( expmat( R * dh ), 1., 1); // normalise row(3rd arg=1) to 1  
		
		// update P & A 
		P  = normalise(clamp(Q.t() * P,0.,1.), 1., 0); // norm cols to 1 
	}
	
	void update_edge_length_substPerSite(const double dh){
		int k; 
		for (int iu = 0; iu < extant.size(); iu++){
			if (extant.at(iu)){
				for (int k = 0; k < m; k++){
					edge_length_substPerSite_nodeWise(iu) += R::rpois( dh * substitutionRates(k) *sequenceLength ); 
				}
			}
		}
	}
	
	//~ List cosim()
	void cosim()
	{
		RNGScope scope;
		//~ find index of first time *after* maxSampleTime
		int start_it = 0; 
		for (start_it = 0; start_it < times.size() - 1; start_it++){
			if (times(start_it) < maxSampleTime){
				start_it--; 
				break;
			}
		}
		// note times decreasing
		for(it = start_it; it < times.size()-1; it++)
		{
			A = sum(P,1); // row sum
			
			t =  times(it);
			t1 =times(it+1);  //t - deltat;
			h = maxSampleTime - t; 
			h0 = h; 
			h1 = maxSampleTime - t1; 
			F = as<mat>(Fs[it]);
			G = as<mat>(Gs[it]);
			Y = clamp( as<vec>(Ys[it]), 1e-12, INFINITY);
			update_states( h1-h0); //
			if (sequenceLength > 0 ){
				update_edge_length_substPerSite( h1 - h0 );
			}
			// co rates
			double lambda = 0.; 
			lambda_uv = zeros(n + Nnode, n + Nnode ); 
			for (i = 0; i < (n + internalNodesAdded); i++){
				for (j = (i+1); j < (n + internalNodesAdded); j++){
					if ( (extant.at(i)) && (extant.at(j)) ){
						lambda_uv(i,j) = co_rate_uv( i, j ); //
						lambda += lambda_uv(i,j);
					}
				}
			}
			int nco = R::rpois( (h1-h0) *  lambda ); 
			// co events
			for (int ico = 0; ico < nco; ico++){
				// ensure that at least two lines are extant if doing co: 
					// find lines to co
					r = R::runif( 0., lambda);
					double rr = 0.; 
					bool founduv = false;
					for (i = 0; i < (n + internalNodesAdded); i++){
						for (j = (i+1); j < (n + internalNodesAdded); j++){
							if ( extant.at(i) && extant.at(j) ){
								rr += lambda_uv(i,j);
								if (rr > r){
									u = i; 
									v = j; 
									founduv=true; 
									break;
								}
							}
						}
						if (rr > r ) break; 
					}
					if (!founduv) {
						break;
					}
					lambda_uv(u, v) = 0.; 
					lambda_uv(v,u) = 0.; 
					
					// do the co
					a = n + internalNodesAdded;
					internalNodesAdded++;
					heights(a) = h;
					
					edge(edgesAdded, 0) = a;
					edge(edgesAdded, 1) = u;
					edge_length(edgesAdded) = heights(a)-heights(u);
					edge_length_substPerSite(edgesAdded) = edge_length_substPerSite_nodeWise(u);
					edgesAdded++;
					edge(edgesAdded, 0) = a;
					edge(edgesAdded, 1) = v;
					edge_length(edgesAdded) = heights(a)-heights(v);
					edge_length_substPerSite(edgesAdded) = edge_length_substPerSite_nodeWise(v);
					edgesAdded++;
					
					//lstate, mstate of a;
					// this gives overly ladder-like trees?
					//~ lstates(a,donordeme) = 1.0;
					//~ mstates(a,donordeme) = 1.0;
					double slsa = 0.;
					lstates.col(a) = zeros(m); 
					for ( int w = 0; w < m ; w++){
						for (int z = 0; z < m; z++){
							lstates(w,a) += F(w,z) * ( P(w,u)*P(z,v) / Y(z)/Y(w) +  P(z,u)*P(w,v) / Y(z)/Y(w) ) ; 
						}
						slsa += lstates(w,a);
					}
					for (int w = 0; w < m; w++){
						lstates(w,a) /= slsa; 
					}
					lstates.col(a) = lstates.col(a) / sum(lstates.col(a)); 
					P.col(a) = lstates.col(a);
					
					extant.at(u) = false;
					extant.at(v) = false;
					extant.at(a) = true;
					numberExtant--;
					P.col(u) = zeros(m);
					P.col(v) = zeros(m);
					
					// set ustates for u and v
					ustates.col(u) = P.col(u);
					ustates.col(v) = P.col(v); 
					
					//update mstates of lines not in co 
					if (finiteSizeCorrection){ //TODO 
						//~ std::cout << " finiteSizeCorrection not implemented " << endl; 
						//~ throw 1;
						finite_size_correction2(a );
					}
					
			}
			// terminate early if all nodes added 
			if ( internalNodesAdded >=(n-1) ) break; 
			// new samples
			incorporate_samples(h0
			  , h1
			);
		}
		//add polytomous root if not all internal nodes added
		if (internalNodesAdded < Nnode) {
			a = n + internalNodesAdded;
			heights(a) = h;
			//~ for (int ico = 0; ico < numberExtant; ico++ ){
			int sumextant= 0;
			for (int ico = 0; ico < (n + Nnode); ico++){
				if (extant.at(ico)){
					sumextant++;
				}
			}
			u = 0; 
			for (int ico = 0; ico < sumextant; ico++ ){
				//find u & v
				bool foundu =  false;
				for (u = u; u < (n + internalNodesAdded); u++){
					if (extant.at(u)){
						extant.at(u) = false; 
						foundu = true; 
						break;
					}
				}
				if (foundu){
					edge(edgesAdded, 0) = a;
					edge(edgesAdded, 1) = u;
					edge_length(edgesAdded) = heights(a)-heights(u);
					edgesAdded++;
				}
			}
		}
		if (samplesAdded < sortedSampleHeights.size()){
			std::cout << "Error: Time axis did not cover all sample times." << std::endl; 
			throw 1; 
		}
		
		if (sequenceLength > 0 ){
			edge_length_substPerSite /= (double)sequenceLength; 
		}
	}

};




//[[Rcpp::export]]
List simulateTreeCpp3x0(const vec times
	  , const List Fs, const List Gs, const List Ys
	  , const vec sortedSampleHeights
	  , const mat sortedSampleStates // m X n 
	  , double maxSampleTime 
	  , const int m
	  , bool finiteSizeCorrection
	  , vec substitutionRates
	  , int sequenceLength)
{
	CoSim cosim(times, Fs, Gs, Ys, sortedSampleHeights, sortedSampleStates, maxSampleTime, m, finiteSizeCorrection
	  ,  substitutionRates
	  ,  sequenceLength);
	//~ return cosim.cosim(); 
	cosim.cosim(); 
	List ret;
	ret["edge"] = cosim.edge;
	ret["edge.length"] = cosim.edge_length;
	ret["edge.length.sps"] = cosim.edge_length_substPerSite; // 
	ret["n"] = cosim.n;
	ret["Nnode"] = cosim.internalNodesAdded;//Nnode
	ret["lstates"] = cosim.lstates;
	ret["ustates"] = cosim.ustates;
	ret["heights"] = cosim.heights;
	ret["samplesAdded"] = cosim.samplesAdded; 
	return ret ;
}

