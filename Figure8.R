#Figure 8a
library( "networkD3" )


Sankey <- read.table( "Sankey.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
Sankey$SM_category <- unlist( lapply( as.character( Sankey$SM_category ), function(x) { return( ifelse( is.na( x ), "Other1", x ) ); } ) ); 
Sankey$TDL_category <- unlist( lapply( as.character( Sankey$TDL_category ), function(x) { return( ifelse( is.na( x ), "Other2", ifelse( x == "", "Other2", x ) ) ); } ) ); 
Sankey$DRUG <- unlist( lapply( as.character( Sankey$DRUG ), function(x) { return( ifelse( is.na( x ), "0", "1" ) ); } ) ); 

print( table( Sankey$SM_category ) )
print( table( Sankey$TDL_category ) )
print( table( Sankey$DRUG ) )

# define "nodes"
dfNode <- data.frame( label=c( "Clinical_Precedence", "Discovery_Precedence", "Predicted_Tractable", "Unknown", "Other1", "Tclin", "Tchem", "Tbio", "Tdark", "Other2", "1", "0" ) ); 
dfNode$name <- as.character( seq( 0, nrow( dfNode ) - 1 ) ); 
dfNode$group <- as.character( dfNode$label ); 

# define "links"
t1 <- as.data.frame( table( Sankey[ , c( "SM_category", "TDL_category" ) ] ) ); 
t2 <- as.data.frame( table( Sankey[ , c( "TDL_category", "DRUG" ) ] ) ); 
colnames( t1 ) <- c( "source", "target", "value" ); 
colnames( t2 ) <- c( "source", "target", "value" ); 
dfLink <- rbind( t1, t2 ); 
dfLink$group <- unlist( apply( dfLink, 1, function(x) { 
	
	if( length( intersect( c( "Discovery_Precedence", "Clinical_Precedence", "Tclin", "1" ), x ) ) > 0 ) { 
		return( "Tclin" ); 
	} else if( length( intersect( c( "Predicted_Tractable", "Tchem" ), x ) ) > 0 ) {
		return( "Tchem" ); 
	} else if( length( intersect( c( "Tbio" ), x ) ) > 0 ) {
		return( "Tbio" ); 
	} else {
		return( "Tdark" ); 
	} 
	
	} ) )
dfLink$source <- as.numeric( as.character( dfNode[ match( as.character( dfLink$source ), as.character( dfNode$label ) ), "name" ] ) ); 
dfLink$target <- as.numeric( as.character( dfNode[ match( as.character( dfLink$target ), as.character( dfNode$label ) ), "name" ] ) ); 
dfLink <- dfLink[ with( dfLink, order( `source`, target ) ), ]; 

my_color <- 'd3.scaleOrdinal().domain( ["Discovery_Precedence", "Clinical_Precedence", "Unknown", "Predicted_Tractable", "Other1", "Tclin", "Tchem", "Tbio", "Tdark", "Other2", "0", "1"] ).range( ["#00acd2", "#ffecd8", "#BFBFBF", "#c7aed6", "#E5E5E5", "#00acd2", "#ffecd8", "#c7aed6", "#BFBFBF", "#E5E5E5", "#BFBFBF", "#00acd2"] )'; 
sankeyNetwork( Links=dfLink, Nodes=dfNode, Source="source", 
	Target="target", Value="value", NodeID="name", NodeGroup="group", LinkGroup="group", colourScale=my_color,
	units="TWh", fontSize=12, nodeWidth=30, height=500, width=500, iterations=0 ); 




