
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
		 , finiteSizeCorrection) {
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
		 , finiteSizeCorrection)
}

colik2cpp <- function(heights, Fs, Gs, Ys, eventIndicator, eventIndicatorNode, eventHeights, sortedSampleStates, daughters, n, Nnode, m, AgtYboundaryCondition) {
	.Call(  'sourceCpp_0_colik2cpp', PACKAGE='phydynR'
	 , heights, Fs, Gs, Ys, eventIndicator, eventIndicatorNode, eventHeights, sortedSampleStates, daughters, n, Nnode, m, AgtYboundaryCondition
	)
}



updateWCpp <- function( W
  , psi_a 
  , utips
  , vtips
  , utipsW
  , vtipsW )
{
	.Call( 'sourceCpp_0_updateWCpp', PACKAGE='phydynR'
	  , W
	  , psi_a 
	  , utips
	  , vtips
	  , utipsW
	  , vtipsW
	)
}


