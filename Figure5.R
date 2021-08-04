#Figure 5a 5b
library( "reshape2" )
library( "ggplot2" )
library( "grid" )

GScore <- read.table( "Gscore_Del.txt", sep="\t", quote="", header=FALSE, check.names=FALSE )
colnames( GScore ) <- unlist( apply( GScore, 2, function(x) {paste(unlist(strsplit(x[1],"_"))[1],x[2],sep="_")} ) )
colnames( GScore )[ 1 ] <- "Gene"
GScore <- GScore[ 3 : nrow( GScore ), ]
GScore.Gene.lst <- GScore[ , "Gene" ]
GScore.Gene.shown.lst <- make.names( GScore.Gene.lst, unique=TRUE )
GScore[ , "Gene" ] <- GScore.Gene.shown.lst
GScore.study.lst <- unique( unlist( lapply( colnames( GScore )[ 2 : ncol( GScore ) ], function(x) {unlist(strsplit(x,"_"))[1]} ) ) )

GScore <- melt( GScore, id="Gene" )
colnames( GScore )[ 2:3 ] <- c( "type", "GScore" )
GScore[ , "study" ] <- unlist( lapply( GScore[ , "type" ], function(x) {unlist(strsplit(as.character(x),"_"))[1]} ) )
GScore[ , "type" ] <- unlist( lapply( GScore[ , "type" ], function(x) {unlist(strsplit(as.character(x),"_"))[2]} ) )

# use ggplot

GScore[ , "type" ] <- unlist( apply( GScore, 1, function(x) {paste(x["study"],x["type"],sep="_")} ) )
study.order <- c( "KIRC","KIRP","KICH","PRAD","BLCA","THYM","MESO","LUAD","LUSC","HNSC","ESCA","STAD","COAD","READ","PAAD","LIHC","CHOL","ACC","PCPG","THCA","BRCA","OV","CESC","UCEC","UCS","SARC","TGCT","LAML","DLBC","LGG","GBM","UVM","SKCM" )
study.order <- unlist( lapply( study.order, function(x) {return( paste( x, c( "Amp", "Del" ), sep="_" ) ); } ) ); 

svg( "TCGA_Del.svg", width=14 * 0.7, height=28 * 0.7 )
p <- ggplot( GScore, aes( x=factor( type, levels=study.order ), y=factor( Gene, levels=rev( GScore.Gene.shown.lst ) ) ) ) + 
	geom_point( aes( size=sqrt( as.numeric( GScore ) ), colour=ifelse( grepl( "Amp", as.character( type ) ), "Amp", "Del" ) ), alpha=0.5, shape=16 ) + 
	xlab( "" ) + 
	ylab( "" ) +
	theme( axis.ticks.length = unit( 0, "cm" ) ) +
	theme( plot.background = element_rect( colour = "white" ), plot.margin = unit( c( 0, 0, 0, 0 ), "mm" ) ) +
	theme( panel.background = element_rect( fill = "white" ), panel.margin = unit( c( 0, 0, 0, 0 ), "mm" ), panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
	theme( strip.background = element_rect( colour = "white" ) ) +
	scale_colour_manual( name="type", values=c( "Amp"="red","Del"="blue" ) ) + 
	# theme_minimal() +
	theme( axis.text.x = element_text( angle=90, hjust=1 ) ) +
	scale_x_discrete( breaks=paste( GScore.study.lst, "Amp", sep="_" ), labels=GScore.study.lst, expand=c( 0.1, 0.1 ) ) +
	scale_size( "GScore", breaks=sqrt( seq( 0.2, 0.9, 0.2 ) ) , labels=seq( 0.2, 0.9, 0.2 ) )
print( p )	
dev.off()

GScore <- read.table( "Gscore_Amp.txt", sep="\t", quote="", header=FALSE, check.names=FALSE )
colnames( GScore ) <- unlist( apply( GScore, 2, function(x) {paste(unlist(strsplit(x[1],"_"))[1],x[2],sep="_")} ) )
colnames( GScore )[ 1 ] <- "Gene"
GScore <- GScore[ 3 : nrow( GScore ), ]
GScore.Gene.lst <- GScore[ , "Gene" ]
GScore.Gene.shown.lst <- make.names( GScore.Gene.lst, unique=TRUE )
GScore[ , "Gene" ] <- GScore.Gene.shown.lst
GScore.study.lst <- unique( unlist( lapply( colnames( GScore )[ 2 : ncol( GScore ) ], function(x) {unlist(strsplit(x,"_"))[1]} ) ) )

GScore <- melt( GScore, id="Gene" )
colnames( GScore )[ 2:3 ] <- c( "type", "GScore" )
GScore[ , "study" ] <- unlist( lapply( GScore[ , "type" ], function(x) {unlist(strsplit(as.character(x),"_"))[1]} ) )
GScore[ , "type" ] <- unlist( lapply( GScore[ , "type" ], function(x) {unlist(strsplit(as.character(x),"_"))[2]} ) )

# use ggplot

GScore[ , "type" ] <- unlist( apply( GScore, 1, function(x) {paste(x["study"],x["type"],sep="_")} ) )
study.order <- c( "KIRC","KIRP","KICH","PRAD","BLCA","THYM","MESO","LUAD","LUSC","HNSC","ESCA","STAD","COAD","READ","PAAD","LIHC","CHOL","ACC","PCPG","THCA","BRCA","OV","CESC","UCEC","UCS","SARC","TGCT","LAML","DLBC","LGG","GBM","UVM","SKCM" )
study.order <- unlist( lapply( study.order, function(x) {return( paste( x, c( "Amp", "Del" ), sep="_" ) ); } ) ); 

svg( "TCGA_Amp.svg", width=14 * 0.7, height=28 * 0.7 )
p <- ggplot( GScore, aes( x=factor( type, levels=study.order ), y=factor( Gene, levels=rev( GScore.Gene.shown.lst ) ) ) ) + 
	geom_point( aes( size=sqrt( as.numeric( GScore ) ), colour=ifelse( grepl( "Amp", as.character( type ) ), "Amp", "Del" ) ), alpha=0.5, shape=16 ) + 
	xlab( "" ) + 
	ylab( "" ) +
	theme( axis.ticks.length = unit( 0, "cm" ) ) +
	theme( plot.background = element_rect( colour = "white" ), plot.margin = unit( c( 0, 0, 0, 0 ), "mm" ) ) +
	theme( panel.background = element_rect( fill = "white" ), panel.margin = unit( c( 0, 0, 0, 0 ), "mm" ), panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) + 
	theme( strip.background = element_rect( colour = "white" ) ) +
	scale_colour_manual( name="type", values=c( "Amp"="red","Del"="blue" ) ) + 
	# theme_minimal() +
	theme( axis.text.x = element_text( angle=90, hjust=1 ) ) +
	scale_x_discrete( breaks=paste( GScore.study.lst, "Amp", sep="_" ), labels=GScore.study.lst, expand=c( 0.1, 0.1 ) ) +
	scale_size( "GScore", breaks=sqrt( seq( 0.2, 0.9, 0.2 ) ) , labels=seq( 0.2, 0.9, 0.2 ) )
print( p )	
dev.off()


#Figure 5c
Percentage <- read.table( "Mut%.txt", sep="\t", quote="", header=TRUE, check.names=FALSE )
Vote <- read.table( "Mut_5way.txt", sep="\t", quote="", header=TRUE, check.names=FALSE )

gene.lst <- as.character( Percentage[ , "gene" ] )
gene.shown.lst <- make.names( gene.lst, unique=TRUE )
Percentage <- cbind( data.frame( gene=gene.shown.lst ), Percentage[ , 3 : ncol( Percentage ) ] )
Vote <- cbind( data.frame( gene=gene.shown.lst ), Vote[ , 3 : ncol( Vote ) ] )

study.lst <- colnames( Percentage )[ 2 : ncol( Percentage ) ]

df <- cbind( melt( Percentage, id="gene" ), melt( Vote, id="gene" )[ , "value", drop=FALSE ] )
colnames( df ) <- c( "gene", "study", "percentage", "vote" )

# use symbol

df[ , "x" ] <- match( df[ , "study" ], study.lst )
df[ , "y" ] <- match( df[ , "gene" ], rev( gene.shown.lst ) )
df[ , "color" ] <- ifelse( df$vote==0, "gainsboro", "darkorange2" )
df[ , "alpha" ] <- ifelse( df$vote==0, 1, df$vote / 5 )
df <- df[ order( df$percentage ), ]
df <- subset( df, percentage !=0 )
df <- subset( df, vote >= 2 )
df[ , "alpha" ] <- (( df$vote ) / 5 )

svg( "TCGA.Mutation.symbol.svg", width=18, height=28 )
layout( matrix( c( 1, 1, 1, 2, 3, 4 ), 3, 2, byrow=FALSE ) )

# par( pin=c( length( as.character( unique( study.lst ) ) ) / 2, length( as.character( unique( gene.shown.lst ) ) ) / 2 ) )
symbols( x=df$x, y=df$y, circles=sqrt( df$percentage ), inches=0.3 / 2, ann=F, bg=alpha( df$color, df$alpha ), fg="white", axes=FALSE )
axis( 1, at=c( 1:length( study.lst ) ), labels=study.lst, tick=FALSE, pos=0, hadj=1, las=3 ) # size = 1 means below
axis( 2, at=c( 1:length( gene.shown.lst ) ), labels=rev( gene.shown.lst ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 

size_legend <- signif( quantile( as.numeric( df$percentage ) ), 1 )
par( pin=c( 7 / 2, 1 / 2 ) )
symbols( x=factor( seq( 1, 7 ) ), y=rep( 1, 7 ), circles=sqrt( c( size_legend[ 1 ], 5, 10, 20, 30, 40, size_legend[ 5 ] ) ), inches=0.3 / 2, ann=F, bg="darkorange2", fg="white", axes=FALSE )
axis( 2, at=1, labels="percentage", tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 
axis( 1, at=c( 1:7 ), labels=paste( as.character( c( size_legend[ 1 ], 5, 10, 20, 30, 40, size_legend[ 5 ] ) ), "%", sep="" ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 1 means below

vote_legend <- c( 1, 2, 3, 4, 5)
par( pin=c( 4 / 2, 1 / 2 ) )
symbols( x=factor( seq( 1, 5 ) ), y=rep( 1, 5 ), circles=rep( 15, 5 ), inches=0.25 / 2, ann=F, bg=alpha( rep( "darkorange2", 5 ), ( vote_legend ) / 5 ), fg="white", axes=FALSE )
axis( 2, at=1, labels="vote", tick=FALSE, pos=0, hadj=1, las=1 ) # vote = 2 means left 
axis( 1, at=c( 1:5 ), labels=as.character( vote_legend ), tick=FALSE, pos=0, hadj=1, las=1 ) # vote = 1 means below

dev.off()


#Figure 5d
Percentage <- read.table( "FusionEvents.txt", sep="\t", quote="", header=TRUE, check.names=FALSE )
Vote <- read.table( "FusionEvents.Vote.txt", sep="\t", quote="", header=TRUE, check.names=FALSE )

gene.lst <- as.character( Percentage[ , "gene" ] )
gene.shown.lst <- make.names( gene.lst, unique=TRUE )
Percentage <- cbind( data.frame( gene=gene.shown.lst ), Percentage[ , 3 : ncol( Percentage ) ] )
Vote <- cbind( data.frame( gene=gene.shown.lst ), Vote[ , 3 : ncol( Vote ) ] )

study.lst <- colnames( Percentage )[ 2 : ncol( Percentage ) ]

df <- cbind( melt( Percentage, id="gene" ), melt( Vote, id="gene" )[ , "value", drop=FALSE ] )
colnames( df ) <- c( "gene", "study", "percentage", "vote" )


# use symbol

df[ , "x" ] <- match( df[ , "study" ], study.lst )
df[ , "y" ] <- match( df[ , "gene" ], rev( gene.shown.lst ) )
df[ , "color" ] <- ifelse( df$vote==1, "gainsboro", "green" )
df[ , "alpha" ] <- ifelse( df$vote==0, 1, df$vote / 5 )
df <- df[ order( df$percentage ), ]
df <- subset( df, percentage >=5 )
df <- subset( df, vote >= 1 )
df[ , "alpha" ] <- ifelse( df$vote < 2, 1, ( df$vote - 1 ) / 4 )

svg( "TCGA.FusionEvents5.symbol.svg", width=20, height=28 )
layout( matrix( c( 1, 1, 1, 2, 3, 4 ), 3, 2, byrow=FALSE ) )

# par( pin=c( length( as.character( unique( study.lst ) ) ) / 2, length( as.character( unique( gene.shown.lst ) ) ) / 2 ) )
symbols( x=df$x, y=df$y, circles=sqrt( df$percentage ), inches=0.4 / 2, ann=F, bg=alpha( df$color, df$alpha ), fg="white", axes=FALSE )
axis( 1, at=c( 1:length( study.lst ) ), labels=study.lst, tick=FALSE, pos=0, hadj=1, las=3 ) # size = 1 means below
axis( 2, at=c( 1:length( gene.shown.lst ) ), labels=rev( gene.shown.lst ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 

size_legend <- signif( quantile( as.numeric( df$percentage ) ), 1 )
# par( pin=c( 5 / 2, 1 / 2 ) )
# symbols( x=factor( seq( 1, 5 ) ), y=rep( 1, 5 ), circles=sqrt( c( 1, 3, 6, 9, 12 ) ), inches=0.4 / 2, ann=F, bg="green", fg="white", axes=FALSE )
# symbols( x=factor( seq( 1, 5 ) ), y=rep( 1, 5 ), circles=sqrt( size_legend ), inches=0.4 / 2, ann=F, bg="green", fg="white", axes=FALSE )
# symbols( x=factor( seq( 1, 5 ) ), y=rep( 1, 5 ), circles=sqrt( c( 0.5, 3, 6, 9, 40 ) ), inches=0.4 / 2, ann=F, bg="green", fg="white", axes=FALSE )
# axis( 2, at=1, labels="percentage", tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 
# axis( 1, at=c( 1:5 ), labels=paste( as.character( c( 0.5, 3, 6, 9, 40 ) ), "%", sep="" ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 1 means below

par( pin=c( 7 / 2, 1 / 2 ) )
symbols( x=factor( seq( 1, 7 ) ), y=rep( 1, 7 ), circles=sqrt( c( 1, 5, 10, 20, 50, 100, size_legend[ 5 ] ) ), inches=0.4 / 2, ann=F, bg="green", fg="white", axes=FALSE )
axis( 2, at=1, labels="percentage", tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 
axis( 1, at=c( 1:7 ), labels=paste( as.character( c(1, 5, 15, 30, 50, 100, size_legend[ 5 ] ) ), "%", sep="" ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 1 means below

color_legend <- c( "gainsboro", "green" )
par( pin=c( length( color_legend ) / 2, 1 / 2 ) )
symbols( x=factor( seq( 1, length( color_legend ) ) ), y=rep( 1, length( color_legend ) ), circles=rep( median( size_legend ), length( color_legend ) ), inches=0.4 / 2, ann=F, bg=color_legend, fg="white", axes=FALSE );
axis( 2, at=1, labels="Color", tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 
axis( 1, at=c( 1:length( color_legend ) ), labels=c( "NS", "S" ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 1 means below

vote_legend <- c( 2, 3, 4, 5 )
par( pin=c( 4 / 2, 1 / 2 ) )
symbols( x=factor( seq( 2, 5 ) ), y=rep( 1, 4 ), circles=rep( 15, 4 ), inches=0.4 / 2, ann=F, bg=alpha( rep( "green", 4 ), ( vote_legend - 1 ) / 4 ), fg="white", axes=FALSE )
axis( 2, at=1, labels="vote", tick=FALSE, pos=0, hadj=1, las=1 ) # vote = 2 means left 
axis( 1, at=c( 1:4 ), labels=as.character( vote_legend ), tick=FALSE, pos=0, hadj=1, las=1 ) # vote = 1 means below

dev.off()






