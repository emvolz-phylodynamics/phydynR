/* 
 * TODO Maybe duplicate likelihood surv terms if co events happen concurrently
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

static const int SAMPLE  = 0; 
static const int CO = 1; 


class CoLik{
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
	
	vec heights;
	
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
	
	// co stuff and tree stuff
	uvec eventIndicator; //(int) SAMPLE or CO
	uvec eventIndicatorNode; //(int) node involved at each event (R index!)
	vec eventHeights; //(double)
	umat daughters; // daughters of each node ((Nnode+n) X 2)
		
	//lik stuff
	double total_lambda; 
	double loglik; 
	vec survTerms; 
	vec rateTerms; 
	
	double AgtYboundaryCondition; // penalises case where A > Y TODO 
	  
	CoLik( vec in_times, List in_Fs,  List in_Gs,  List in_Ys //NOTE times decreasing
	  , vec in_sortedSampleHeights
	  , mat in_sortedSampleStates // m X n 
	  , uvec in_eventIndicator // sample or co (int)
	  , uvec in_eventIndicatorNode // node involved at each event
	  , vec in_eventHeights
	  , umat in_daughters // daughters of each node
	  , double in_maxSampleTime 
	  , int in_m
	  , bool in_finiteSizeCorrection
	  , double in_maxHeight) : times(in_times),Fs(in_Fs),Gs(in_Gs),Ys(in_Ys),sortedSampleHeights(in_sortedSampleHeights),sortedSampleStates(in_sortedSampleStates),maxSampleTime(in_maxSampleTime),m(in_m),finiteSizeCorrection(in_finiteSizeCorrection)
	  ,eventIndicator(in_eventIndicator)
	  ,eventIndicatorNode(in_eventIndicatorNode)
	  ,eventHeights(in_eventHeights)
	  ,daughters(in_daughters)
	  ,maxHeight(in_maxHeight)
	{
		
		n = sortedSampleHeights.size();
		Nnode = sum(eventIndicator);//n -1; 
//~ cout << "Nnode " << Nnode << endl;
		// note times decreasing !
		tfin = times(times.size()-1);
		maxHeight = maxSampleTime -tfin;
		
		//~ extant(n + Nnode, false); 
		extant.resize(n + Nnode); std::fill( extant.begin(), extant.end(), false); 
		numberExtant = 0;
		
		//lik stuff
		total_lambda = 0.;
		loglik = 0.;
		survTerms.resize(times.size()); std::fill( survTerms.begin(), survTerms.end(), 0.);
		rateTerms.resize(n+Nnode); std::fill( rateTerms.begin(), rateTerms.end(), 0.); 
		
		heights = zeros(n+Nnode); heights.fill(INFINITY);
		
		lambda_uv = zeros(n + Nnode, n + Nnode ); 
		P = zeros( m, n + Nnode );
		lstates = zeros(m, n + Nnode);
		ustates = zeros(m, n + Nnode);
		
		internalNodesAdded = 0;
		samplesAdded = 0; // counter for terminal branches
		
		A = zeros(m); 
		F = zeros(m,m);
		G = zeros(m,m);
		
		h = 0.; 
	}
	
	void incorporate_samples(const int a)
	{
		int i = a - 1; 
		extant.at(i) = true; 
		heights(i)  =sortedSampleHeights(i);
		P.col(i) = sortedSampleStates.col(i); 
		lstates.col(i) = sortedSampleStates.col(i);
		samplesAdded++;
		numberExtant++;
	}
	
	
	void finite_size_correction2(int a)
	{ //
		int ia = a - 1;
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
			if (extant.at(iu) && (iu != ia)){
				p_u = P.col(iu); 
				//~ Amin_pu = clamp(( A - p_u), 1., INFINITY ); 
				Amin_pu = clamp(( A - p_u), 1e-3, INFINITY ); 
				rterm = P.col(ia) / Amin_pu ; 
				rho = A / Amin_pu; 
				lterm = dot( rho, P.col(ia)); //
				p_u = p_u % clamp((lterm - rterm), 0., INFINITY) ; // l > r
				p_u = p_u / sum(p_u ) ; 
				P.col(iu) = p_u; 
			}
		}
	}	
	
	void update_states( const double dh){
		vec A_Y = clamp(A / Y, 0., 1.); 
		
		// make R & Q 
		//~ (F(k,l) + G(k,l)) *  Q(x, l,z)/  std::max(Q(x,l,z), Y(l));
		mat R = (F+G).t();
		//~ R  = R.each_col() / Y; 
		R.each_col() /= Y; 
		R.diag().zeros();
		R.diag() = -sum(R,1) - ((F *  A_Y) / Y) ;
		mat Q = normalise( expmat( R * dh ), 1., 1); // normalise row(3rd arg=1) to 1  
		
		// update P & A 
		P  = normalise(clamp(Q.t() * P,0.,1.), 1., 0); // norm cols to 1 
	}
	
	double co_rate_uv(const int i, const int j)
	{ //TODO unit test
		double lambdauv = 0.; 
		for (k = 0; k < m; k++){
			for (l = 0; l < m; l++){
				lambdauv += F(k,l)*( P(k,i)*P(l,j) + P(k,j)*P(l,i))/( Y(k)*Y(l) );
			}
		}
		return lambdauv; 
		//~ vec pu__Y = clamp(P.col(i),0., 1.)/clamp(Y, 1e-12, INFINITY) ;
		//~ vec pv__Y = clamp(P.col(j), 0., 1.)/clamp(Y, 1e-12, INFINITY);
		//~ return as_scalar(dot(pu__Y , (F * pv__Y)) + dot(pv__Y , (F * pu__Y))); 
	}
	
	void update_lambda_uv(){
		//~ lambda_uv = zeros(n + Nnode, n + Nnode ); 
		total_lambda = 0.; 
		A = sum(P,1); // row sum
		for (i = 0; i < (n + internalNodesAdded); i++){
			for (j = (i+1); j < (n + internalNodesAdded); j++){
				if ( (extant.at(i)) && (extant.at(j)) ){
					lambda_uv(i,j) = co_rate_uv( i, j ); //
					total_lambda += lambda_uv(i,j);
				}
			}
		}
	}
	
	
	void colik()
	{
		RNGScope scope;
		int eventsPassed = 0;
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
			//~ A = sum(P,1); // row sum
			t =  times(it);
			t1 =times(it+1);  //t - deltat;
			h = maxSampleTime - t; 
			h0 = h; 
			h1 = maxSampleTime - t1; 
			
			if (h0 >=maxHeight){
				break;
			}
			
			F = as<mat>(Fs[it]);
			G = as<mat>(Gs[it]);
			Y = clamp( as<vec>(Ys[it]), 1e-12, INFINITY);
//~ cout << h0 << endl; 
//~ cout << F << endl; 
//~ cout << Y << endl; 
//~ cout  << endl; 

			update_states( h1-h0); //
//~ cout << " update states "<< endl; 
			// co rates
			//~ double lambda = 0.; 
			update_lambda_uv();
//~ cout << " update lambdauv " << endl ;
			// events in (h0,h1)
			double hevent;
			int eventType;
			for (int ievent = eventsPassed; ievent < eventIndicator.size(); ievent++){
				hevent = eventHeights(ievent);
				eventType = eventIndicator(ievent);
				a = eventIndicatorNode(ievent); 
//~ cout << ievent << " "<<hevent << " " << eventType << " " << a << " " << h0 << " " << h1 << endl ;
				if (hevent > h1) {
					break;
				}
				if (hevent > h0 && hevent <= h1){
					if (eventType == SAMPLE){
//~ cout << "inc samp " << a << endl; 
						incorporate_samples(a);
//~ cout << "done inc samp " << a << endl; 
						update_lambda_uv(); //TODO make a faster version of this just for a
//~ cout << "updated lambdauv " << total_lambda << endl; 
					} else if (eventType==CO){
//~ cout << " co event " << endl; 
						u = daughters(a-1,0); 
						v = daughters(a-1,1);
if (!extant.at(u-1) || !extant.at(v-1)){
	cout << "u or v not extant" << endl;
	throw 99; 
}
						double slsa = 0.;
						lstates.col(a-1) = zeros(m); 
						for ( int w = 0; w < m ; w++){
							for (int z = 0; z < m; z++){
								lstates(w,a-1) += F(w,z) * ( P(w,u-1)*P(z,v-1) / Y(z)/Y(w) +  P(z,u-1)*P(w,v-1) / Y(z)/Y(w) ) ; 
							}
							slsa += lstates(w,a-1);
						}
						for (int w = 0; w < m; w++){
							lstates(w,a-1) /= slsa; 
						}
						lstates.col(a-1) = lstates.col(a-1) / sum(lstates.col(a-1)); 
						P.col(a-1) = lstates.col(a-1);
						
						// set ustates for u and v
						ustates.col(u-1) = P.col(u-1);
						ustates.col(v-1) = P.col(v-1);
						
						double lambdauv = co_rate_uv(u-1, v-1); 
						loglik += log( lambdauv ) ;
						rateTerms(a-1) = log( lambdauv ) ;
						//~ loglik += log(lambda_uv(u-1,v-1) ) ;
						//~ rateTerms(a-1) = log(lambda_uv(u-1,v-1) ) ;
						
						extant.at(u-1) = false;
						extant.at(v-1) = false;
						extant.at(a-1) = true;
						numberExtant--;
						P.col(u-1) = zeros(m);
						P.col(v-1) = zeros(m);
						
						internalNodesAdded++;
						
						update_lambda_uv(); //TODO make a faster version of this just for a
						//note ^ updates A
						
						//update mstates of lines not in co 
						if (finiteSizeCorrection){ 
							//~ std::cout << " finiteSizeCorrection not implemented " << endl; 
							//~ throw 1;
							finite_size_correction2(a );
						}
					} 
					eventsPassed++;
				}
			}
			
			loglik -= total_lambda * (h1-h0); 
//~ cout << " loglik " << loglik << endl; 
			survTerms(it) = -total_lambda * (h1-h0); 
//~ cout << " surv " << survTerms(it) << " "<< total_lambda << endl; 
			
			// terminate early if all nodes added 
			if ( internalNodesAdded >= Nnode ){
				break; 
			}
			
		}
		
	}

};




//[[Rcpp::export]]
List colik4cpp(const vec in_times
	  , const List in_Fs
	  , const  List in_Gs
	  , const List in_Ys //NOTE times decreasing
	  , const vec in_sortedSampleHeights
	  , const mat in_sortedSampleStates // m X n 
	  , const uvec in_eventIndicator // sample or co (int)
	  , const uvec in_eventIndicatorNode // node involved at each event
	  , const vec in_eventHeights
	  , const umat in_daughters // daughters of each node
	  , const double in_maxSampleTime 
	  , const int in_m
	  , const bool in_finiteSizeCorrection
	  , const double in_maxHeight
	  )
{
//~ cout << in_sortedSampleHeights << endl;
//~ cout << in_eventIndicator << endl;
//~ cout << in_eventIndicatorNode << endl;
//~ cout << in_daughters << endl;
	CoLik colik(in_times
	  , in_Fs
	  , in_Gs
	  , in_Ys //NOTE times decreasing
	  , in_sortedSampleHeights
	  , in_sortedSampleStates // m X n 
	  , in_eventIndicator // sample or co (int)
	  , in_eventIndicatorNode // node involved at each event
	  , in_eventHeights
	  , in_daughters // daughters of each node
	  , in_maxSampleTime 
	  , in_m
	  , in_finiteSizeCorrection
	  , in_maxHeight);
	 
	colik.colik(); 
	List ret;
	ret["n"] = colik.n;
	ret["Nnode"] = colik.internalNodesAdded;//Nnode
	ret["lstates"] = colik.lstates;
	ret["ustates"] = colik.ustates;
	//~ ret["heights"] = cosim.heights;
	ret["samplesAdded"] = colik.samplesAdded; 
	ret["loglik"] = colik.loglik; 
	ret["rateTerms"] = colik.rateTerms; 
	ret["survTerms"] = colik.survTerms; 
	return ret ;
}

