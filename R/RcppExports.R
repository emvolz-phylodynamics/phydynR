
dAL <- function(t, y, parms) {
	.Call('sourceCpp_0_dAL',  PACKAGE='phydynR', t, y, parms)
}

simulateTreeCpp2 <- function(times
		 ,  Fs
		 ,  Gs
		 ,  Ys
		 , As
		 , sortedCoHeights 
		 , sortedSampleHeights 
		 , sortedSampleStates
		 , maxSampleTime
		 , m
		 , finiteSizeCorrection
		 , DEMES) {
	.Call('sourceCpp_2_simulateTreeCpp2',  PACKAGE='phydynR', times
		 ,  Fs
		 ,  Gs
		 ,  Ys
		 , As
		 , sortedCoHeights 
		 , sortedSampleHeights 
		 , sortedSampleStates
		 , maxSampleTime
		 , m
		 , finiteSizeCorrection
		 , DEMES)
}

simulateTreeCpp3x0 <- function(times
		 ,  Fs
		 ,  Gs
		 ,  Ys
		 , sortedSampleHeights 
		 , sortedSampleStates
		 , maxSampleTime
		 , m
		 , finiteSizeCorrection
		 , substitutionRates
		 , sequenceLength
) {
	.Call('sourceCpp_simulateTreeCpp3x0',  PACKAGE='phydynR'
		 , times
		 ,  Fs
		 ,  Gs
		 ,  Ys
		 , sortedSampleHeights 
		 , sortedSampleStates
		 , maxSampleTime
		 , m
		 , finiteSizeCorrection
		 , substitutionRates
		 , sequenceLength
	)
}



#########################

sourceAttribMultiDemeCpp <- function(heights, Fs, Gs, Ys, eventIndicator, eventIndicatorNode, eventHeights, sortedSampleStates, daughters, n, Nnode, m, AgtYboundaryCondition, maxHeight) {
	.Call(  'sourceCpp_0_sourceAttribMultiDemeCpp', PACKAGE='phydynR'
	 , heights, Fs, Gs, Ys, eventIndicator, eventIndicatorNode, eventHeights, sortedSampleStates, daughters, n, Nnode, m, AgtYboundaryCondition, maxHeight
	)
}


sourceAttribMultiDemeCpp2 <- function(heights, Fs, Gs, Ys, eventIndicator, eventIndicatorNode, eventHeights, sortedSampleStates, daughters, n, Nnode, m, AgtYboundaryCondition, maxHeight, step_size_res) {
	.Call(  'sourceCpp_1_sourceAttribMultiDemeCpp2', PACKAGE='phydynR'
	 , heights, Fs, Gs, Ys, eventIndicator, eventIndicatorNode, eventHeights, sortedSampleStates, daughters, n, Nnode, m, AgtYboundaryCondition, maxHeight, step_size_res
	)
}

#~ sourceCpp_1_sourceAttribMultiDemeCpp2


############# colik modular stuff
rco_finite_size_correction2 <- function(a, p_a, A, extantLines, mstates) {
	.Call( 'sourceCpp_0_rco_finite_size_correction2', PACKAGE='phydynR'
	, a, p_a, A, extantLines, mstates)
}
eventTimes2extant <- function(eventTimes, nodeheights, parentheights) {
	.Call( 'sourceCpp_0_eventTimes2extant' , PACKAGE='phydynR'
	, eventTimes, nodeheights, parentheights)
}
update_alpha0 <- function(pu, pv, Fmat, Yvec, A) {
	.Call(  'sourceCpp_0_update_alpha0', PACKAGE='phydynR'
	,pu, pv, Fmat, Yvec, A)
}
update_states0 <- function(mstates, Q) {
	.Call(  'sourceCpp_0_update_states0', PACKAGE='phydynR'
	 , mstates, Q)
}
update_states1 <- function(mstates, Q, el) {
	.Call(  'sourceCpp_0_update_states1', PACKAGE='phydynR'
	 , mstates, Q, el)
}
solveQALboost0 <- function(times, Fs, Gs, Ys, h0, h1, L0, A0) {
	.Call( 'sourceCpp_2_solveQALboost0' , PACKAGE='phydynR'
	 , times, Fs, Gs, Ys, h0, h1, L0, A0)
}


dQLcpp <- function ( x, FF, GG, YY, A0 ){
	.Call( 'scpp_dQL', PACKAGE='phydynR'
	 , x, FF , GG , YY, A0 )
}

############################

solvePikL0 <- function( times, Fs, Gs, Ys, h0, h1, pik0, step_size_res = 10) {
	.Call( 'sourceCpp_1_solvePikL0', PACKAGE='phydynR', times, Fs, Gs, Ys, h0, h1, pik0, step_size_res )
}


############################

solveQfwd0 <- function(times, Fs, Gs, Ys, deaths, m, h1, h0) {
	.Call('sourceCpp_1_solveQfwd0', PACKAGE='phydynR', times, Fs, Gs, Ys, deaths, m, h1, h0 )
}

