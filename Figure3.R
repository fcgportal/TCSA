#Figure 3d
library( "ggplot2" )
library( "reshape2" )
library( "data.table" )


geneinfo <- read.table( "geneinfo.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 


system.time( expPan <- fread( "TcgaTargetGtex_fpkm", data.table=FALSE ) )
row.names( expPan ) <- as.character( expPan[ , 1 ] ); 
expPan <- expPan[ , 2 : ncol( expPan ) ]; 

phenotype <- read.table( "TcgaTargetGTEX_phenotype.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="" )
colnames( phenotype ) <- gsub( "^_", "", colnames( phenotype ) ); 

# do not use TARGET
phenotype <- subset( phenotype, study %in% c( "GTEX", "TCGA" ) ); 
expPan <- expPan[ , as.character( phenotype$sample ) ]; 

# do not use Testis, Bone Marrow, Lymphocytes or Fibroblasts from GTEx
phenotype <- rbind( subset( phenotype, study == "TCGA" ), subset( phenotype, study == "GTEX" & sample_type != "Cell Line" & primary_site != "Testis" ) ); 
expPan <- expPan[ , as.character( phenotype$sample ) ]; 

# define normal group and each cancer type group
gdc_cases <- read.table( "gdc_cases_update04072019.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="" )
phenotype$Group <- unlist( lapply( as.character( phenotype$sample ), function(x) { 
	
	if( grepl( "^GTEX", x ) ) { 
		
		site <- as.character( subset( phenotype, sample == x )$primary_site ); 
		if( site == "" ) { 
			tissue <- as.character( subset( phenotype, sample == x )$`primary disease or tissue` ); 
			site <- c( "Skin", "Stomach", "Esophagus" )[ match( tissue, c( "Skin - Sun Exposed (Lower Leg)", "Stomach", "Esophagus - Mucosa" ) ) ]; 
		}
		return( paste( "NormalGTEX", site, sep=", " ) ); 
		
	}
	
	t <- unlist( strsplit( x, "-" ) )[ 4 ]; 
	p <- paste( unlist( strsplit( x, "-" ) )[ 1:3 ], collapse="-" ); 
	tt <- as.character( subset( gdc_cases, submitter_id == p )$project_id ); 
	if( grepl( "^1", t ) ) { return( paste( "NormalTCGA", as.character( subset( phenotype, sample == x )$primary_site ), sep=", " ) ); }
	if( grepl( "^2", t ) ) { return( NA ); }
	if( length( tt ) == 0 ) { stop( x, "\n" ); }
	return( tt ); 
	
	} ) )
tab <- table( sort( as.character( phenotype$Group ) ) ); 

row.names( expPan ) <- unlist( lapply( strsplit( row.names( expPan ), "\\." ), function(x) {x[1]} ) )


expMedian <- read.table( "expMedian.txt", header=TRUE, check.names=FALSE, row.names=1, sep="\t", quote="", comment.char="" )

tlst <- colnames( expMedian ); 
tlst0 <- tlst[ grepl( "^NormalGTEX", tlst ) ]
tlst1 <- tlst[ grepl( "^TCGA", tlst ) ]; 

tlst0 <- c( "NormalGTEX, Kidney", "NormalGTEX, Prostate", "NormalGTEX, Bladder", "NormalGTEX, Lung", "NormalGTEX, Esophagus", "NormalGTEX, Stomach", "NormalGTEX, Colon", "NormalGTEX, Small Intestine", "NormalGTEX, Pancreas", "NormalGTEX, Liver", "NormalGTEX, Adrenal Gland", "NormalGTEX, Thyroid", "NormalGTEX, Breast", "NormalGTEX, Fallopian Tube", "NormalGTEX, Ovary", "NormalGTEX, Cervix Uteri", "NormalGTEX, Uterus", "NormalGTEX, Vagina", "NormalGTEX, Brain", "NormalGTEX, Nerve", "NormalGTEX, Blood", "NormalGTEX, Skin", "NormalGTEX, Adipose Tissue", "NormalGTEX, Blood Vessel", "NormalGTEX, Heart", "NormalGTEX, Muscle", "NormalGTEX, Pituitary", "NormalGTEX, Salivary Gland", "NormalGTEX, Spleen" ); 
tlst1 <- c( "TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH", "TCGA-PRAD", "TCGA-BLCA", "TCGA-THYM", "TCGA-MESO", "TCGA-LUAD", "TCGA-LUSC", "TCGA-HNSC", "TCGA-ESCA", "TCGA-STAD", "TCGA-COAD", "TCGA-READ", "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL", "TCGA-ACC", "TCGA-PCPG", "TCGA-THCA", "TCGA-BRCA", "TCGA-OV", "TCGA-CESC", "TCGA-UCEC", "TCGA-UCS", "TCGA-SARC", "TCGA-TGCT", "TCGA-LAML", "TCGA-DLBC", "TCGA-LGG", "TCGA-GBM", "TCGA-UVM", "TCGA-SKCM" ); 

df <- read.table( "caGESPs.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
print( dim( df ) ); # print( showd( df ) ); # 628
print( length( unique( as.character( df$`ENSEMBL Gene ID` ) ) ) ); # 409
print( table( df$`caGESP Tier` ) ); 

df$Tier <- c( "L1", "L2", "L3" )[ match( df$`caGESP Tier`, c( "Tier I", "Tier II", "Tier III" ) ) ]; 
df$Tier <- factor( as.character( df$Tier ), levels=c( "L1", "L2", "L3" ) ); 
print( table( df$`Tier` ) ); 
#     L1     L2     L3
#  170 251 207 

ct_code <- read.table( "Cancertype.code.txt", header=FALSE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
df$CancerType <- paste0( "TCGA-", as.character( ct_code$`V3`[ match( df$`Cancer Type`, ct_code$V2 ) ] ) ); 
print( table( is.na( df$`CancerType` ) ) ); 
print( table( df$`CancerType` ) ); 


# === order SP genes, first according to different criteria, second according to cancer type, then according to SPM
dfP <- df; 
dfP$CancerType <- factor( as.character( dfP$CancerType ), levels=tlst1 ); 
dfP <- dfP[ with( dfP, order( Tier, CancerType ) ), ]; 

spm <- read.table( "spm.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" )
spm$ensg <- unlist( lapply( strsplit( as.character( spm$Acc ), "\\." ), function(x) {x[1]} ) )
row.names( spm ) <- unlist( apply( spm, 1, function(x) { return( paste( x[ "CancerType" ], x[ "ensg" ], sep="_" ) ); } ) )
dfP$identifier <- unlist( apply( dfP, 1, function(x) { return( paste( x[ "CancerType" ], x[ "ENSEMBL Gene ID" ], sep="_" ) ); } ) )
dfP <- cbind( dfP, spm[ as.character( dfP$identifier ) , "SPM", drop=FALSE ] ); 

dfP$ensg <- as.character( dfP$`ENSEMBL Gene ID` ); 
t <- lapply( unique( as.character( dfP$ensg ) ), function(g) { 
	
	t <- subset( dfP, ensg == g ); 
	if( nrow( t ) == 1 ) { return( t ); }
	
	t <- subset( t, Tier == as.character( t[ 1, "Tier" ] ) ); 
	if( nrow( t ) == 1 ) { return( t ); }
	
	t <- t[ with( t, rev( order( SPM ) ) ), ]; 
	t <- subset( t, SPM == t[ 1, "SPM" ] ); 
	if( nrow( t ) == 1 ) { return( t ); }
	
	t <- t[ with( t, order( Tier, CancerType ) ), ]; 
	
	if( is.na( t[ 1, "Tier" ] ) ) { g <<- g; stop( g ); }
	return( t[ 1, ] ); 
	
	} )
dfPnr <- do.call( rbind, t ); 
dfPnr <- dfPnr[ with( dfPnr, order( Tier, CancerType, - SPM ) ), ]; 
dfPnr$yPos <- rev( seq( 1, nrow( dfPnr ) ) ); 
print( table( dfPnr$Tier ) ); 
#  L1  L2  L3 
# 135 161 113 # total is 409


# transform FPKM to fractional density 
dfFrac <- expMedian[ unique( as.character( dfPnr$ensg ) ), c( tlst0, tlst1 ) ]; 
tmp <- apply( dfFrac, 1, function(x) { return( x / sum( x ) ); } ); 
dfFrac <- data.frame( t( tmp ), check.names=FALSE ); 
print( nrow( dfFrac ) )
# [1] 409

dfPlot <- melt( cbind( data.frame( Gene=row.names( dfFrac ) ), dfFrac ), id="Gene" );
dfPlot$yPos <- dfPnr[ match( as.character( dfPlot$Gene ), as.character( dfPnr$ensg ) ), "yPos" ];

dfPlot$variable <- factor( as.character( dfPlot$variable ), levels=c( tlst0, tlst1 ) )
dfPlot2 <- dfPlot; 
dfPlot2$Tier <- as.character( dfPnr[ match( as.character( dfPlot2$Gene ), as.character( dfPnr$ensg ) ), "Tier" ] ); 

p1 <- ggplot() + 
	geom_tile( data=dfPlot, aes( variable, yPos, fill=value ), color=NA ) +
	scale_fill_gradientn( colours=c( "#FFFFFF", "#FFF2F2", "#FFE6E6", "#FFB2B2", "#FF8080", "#FF6666", "#FF4D4D", "#FF3333", "#FF1919", "#FF0000" ), na.value="#E5E5E5", limits=c( 0, 1 ) ) +
	geom_segment( data=data.frame( y=c( unlist( lapply( c( "L1", "L2", "L3" ), function(t) { return( min( subset( dfPnr, Tier == t )$yPos ) - 0.5 ); } ) ), 
		0.5, nrow( dfPnr ) + 0.5 ) ), aes( x=0.5, xend=length( c( tlst0, tlst1 ) ) + 0.5, y=y, yend=y ), colour="#D9D9D9", size=0.4 ) + 
	geom_segment( data=data.frame(), aes( x=length( tlst0 ) + 0.5, xend=length( tlst0 ) + 0.5, y=0.5, yend=nrow( dfPnr ) + 0.5 ), colour="#D9D9D9", size=0.4 ) + 
	scale_x_discrete( breaks=c( tlst0, tlst1 ), labels=gsub( "NormalGTEX, ", "", c( tlst0, tlst1 ) ), position="top" ) +
	scale_y_continuous( breaks=NULL, limits=c( 0.5, nrow( dfPnr ) + 0.5 ), expand=c( 0, 0 ) ) +
	theme( legend.position="right", legend.title=element_blank(), aspect.ratio=1,
		panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_rect( colour="black", fill=NA, size=1.5 ), panel.background=element_blank(),
		axis.line.x=element_blank(), axis.line.y=element_blank(), 
		axis.text.x=element_text( angle=90, hjust=0, vjust=0.5, size=7 ), axis.text.y=element_blank(), 
		axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
		plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
	xlab( "" ) + ylab( "" ) + labs( title=paste( "membrane SP (n=", nrow( dfFrac ), ")", sep="" ) )
ggsave( file="tmp1.svg", plot=p1, width=7, height=7 )
p2 <- ggplot() + 
	geom_tile( data=subset( dfPlot2, Tier == "L1" ), aes( -1, yPos ), fill="#843c0c", color=NA ) +
	geom_tile( data=subset( dfPlot2, Tier == "L2" ), aes( -1, yPos ), fill="#ed7d31", color=NA ) +
	geom_tile( data=subset( dfPlot2, Tier == "L3" ), aes( -1, yPos ), fill="#f4b183", color=NA ) +
	scale_y_continuous( breaks=NULL, limits=c( 0.5, nrow( dfPnr ) + 0.5 ), expand=c( 0, 0 ) ) +
	theme( legend.position="right", legend.title=element_blank(), aspect.ratio=10,
		panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(),
		axis.line.x=element_blank(), axis.line.y=element_blank(), 
		axis.text.x=element_text( angle=90, hjust=0, vjust=0.5, size=7 ), axis.text.y=element_blank(), 
		axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
		plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
	xlab( "" ) + ylab( "" ) + labs( title="" )
ggsave( file="Fig3d.svg", plot=p2, width=3, height=7 )

#Figure 3e
library( "ggplot2" )
library( "reshape2" )
library( "data.table" )


geneinfo <- read.table( "geneinfo.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 


system.time( expPan <- fread( "TcgaTargetGtex_fpkm", data.table=FALSE ) )
row.names( expPan ) <- as.character( expPan[ , 1 ] ); 
expPan <- expPan[ , 2 : ncol( expPan ) ]; 

phenotype <- read.table( "TcgaTargetGTEX_phenotype.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="" )
colnames( phenotype ) <- gsub( "^_", "", colnames( phenotype ) ); 

# do not use TARGET
phenotype <- subset( phenotype, study %in% c( "GTEX", "TCGA" ) ); 
expPan <- expPan[ , as.character( phenotype$sample ) ]; 

# do not use Testis, Bone Marrow, Lymphocytes or Fibroblasts from GTEx
phenotype <- rbind( subset( phenotype, study == "TCGA" ), subset( phenotype, study == "GTEX" & sample_type != "Cell Line" & primary_site != "Testis" ) ); 
expPan <- expPan[ , as.character( phenotype$sample ) ]; 

# define normal group and each cancer type group
gdc_cases <- read.table( "gdc_cases_update04072019.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="" )
phenotype$Group <- unlist( lapply( as.character( phenotype$sample ), function(x) { 
	
	if( grepl( "^GTEX", x ) ) { 
		
		site <- as.character( subset( phenotype, sample == x )$primary_site ); 
		if( site == "" ) { 
			tissue <- as.character( subset( phenotype, sample == x )$`primary disease or tissue` ); 
			site <- c( "Skin", "Stomach", "Esophagus" )[ match( tissue, c( "Skin - Sun Exposed (Lower Leg)", "Stomach", "Esophagus - Mucosa" ) ) ]; 
		}
		return( paste( "NormalGTEX", site, sep=", " ) ); 
		
	}
	
	t <- unlist( strsplit( x, "-" ) )[ 4 ]; 
	p <- paste( unlist( strsplit( x, "-" ) )[ 1:3 ], collapse="-" ); 
	tt <- as.character( subset( gdc_cases, submitter_id == p )$project_id ); 
	if( grepl( "^1", t ) ) { return( paste( "NormalTCGA", as.character( subset( phenotype, sample == x )$primary_site ), sep=", " ) ); }
	if( grepl( "^2", t ) ) { return( NA ); }
	if( length( tt ) == 0 ) { stop( x, "\n" ); }
	return( tt ); 
	
	} ) )
tab <- table( sort( as.character( phenotype$Group ) ) ); 
write.table( as.data.frame( tab ), file="tmp2.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 

row.names( expPan ) <- unlist( lapply( strsplit( row.names( expPan ), "\\." ), function(x) {x[1]} ) )


# load medianFPKM matrix generated by "C:\Users\linzhang\Box Sync\Jiao0_work\2019 Apr 20 Cancer Specificity\Tools\0.MyAnalysis\SP.by5.correct.R"
expMedian <- read.table( "expMedian.txt", header=TRUE, check.names=FALSE, row.names=1, sep="\t", quote="", comment.char="" )

tlst <- colnames( expMedian ); 
tlst0 <- tlst[ grepl( "^NormalGTEX", tlst ) ]
tlst1 <- tlst[ grepl( "^TCGA", tlst ) ]; 
print( length( tlst0 ) ); # 29
print( length( tlst1 ) ); # 33

tlst0 <- c( "NormalGTEX, Kidney", "NormalGTEX, Prostate", "NormalGTEX, Bladder", "NormalGTEX, Lung", "NormalGTEX, Esophagus", "NormalGTEX, Stomach", "NormalGTEX, Colon", "NormalGTEX, Small Intestine", "NormalGTEX, Pancreas", "NormalGTEX, Liver", "NormalGTEX, Adrenal Gland", "NormalGTEX, Thyroid", "NormalGTEX, Breast", "NormalGTEX, Fallopian Tube", "NormalGTEX, Ovary", "NormalGTEX, Cervix Uteri", "NormalGTEX, Uterus", "NormalGTEX, Vagina", "NormalGTEX, Brain", "NormalGTEX, Nerve", "NormalGTEX, Blood", "NormalGTEX, Skin", "NormalGTEX, Adipose Tissue", "NormalGTEX, Blood Vessel", "NormalGTEX, Heart", "NormalGTEX, Muscle", "NormalGTEX, Pituitary", "NormalGTEX, Salivary Gland", "NormalGTEX, Spleen" ); 
tlst1 <- c( "TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH", "TCGA-PRAD", "TCGA-BLCA", "TCGA-THYM", "TCGA-MESO", "TCGA-LUAD", "TCGA-LUSC", "TCGA-HNSC", "TCGA-ESCA", "TCGA-STAD", "TCGA-COAD", "TCGA-READ", "TCGA-PAAD", "TCGA-LIHC", "TCGA-CHOL", "TCGA-ACC", "TCGA-PCPG", "TCGA-THCA", "TCGA-BRCA", "TCGA-OV", "TCGA-CESC", "TCGA-UCEC", "TCGA-UCS", "TCGA-SARC", "TCGA-TGCT", "TCGA-LAML", "TCGA-DLBC", "TCGA-LGG", "TCGA-GBM", "TCGA-UVM", "TCGA-SKCM" ); 


# load SP gene list, by Zhongyi, after filtering out immune-related genes
df <- read.table( "SP.by5.L3.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 

	# do not consider subtype 
dfP <- subset( df, grepl( "^TCGA", CancerType ) );

dfP$Tier <- factor( as.character( dfP$Tier ), levels=c( "L1", "L2", "L3" ) ); 
dfP$CancerType <- factor( as.character( dfP$CancerType ), levels=tlst1 ); 
dfP <- dfP[ with( dfP, order( Tier, CancerType ) ), ]; 


# plot examples 
plotEg0 <- function( GeneA ) { 
	
	# GeneA <- "CLDN6"; 
	print( GeneA ); 
	
	dfPlot <- subset( phenotype, Group %in% c( tlst0, tlst1 ) ); 
	dfPlot$Class <- unlist( lapply( as.character( dfPlot$Group ), function(x) { return( ifelse( grepl( "^TCGA", x ), 1, 0 ) ); } ) );
	
	if( GeneA %in% row.names( expPan ) ) { dfPlot$expA <- as.vector( as.matrix( expPan[ GeneA, as.character( dfPlot$sample ) ] ) ); 
	} else { GeneA <- as.character( subset( geneinfo, gene == GeneA )$ensg ); 
		dfPlot$expA <- as.vector( as.matrix( expPan[ GeneA, as.character( dfPlot$sample ) ] ) ); }
	print( GeneA ); 
	
	xTick <- sort( unique( as.character( dfPlot$Group ) ) ); 
	
	t <- subset( dfP, ensg == GeneA ); 
	if( nrow( t ) > 0 ) { fn <- unlist( apply( t[ 1, ], 1, function(x) { return( paste( "t", x[ "gene" ], x[ "Tier" ], sep=" " ) ); } ) ); 
	} else { fn <- unlist( apply( t[ 1, ], 1, function(x) { return( paste( "t", GeneA, x[ "Tier" ], sep=" " ) ); } ) ); }
	
	ymax <- max( unlist( lapply( c( tlst0, tlst1 ), function(x) { 
		qt <- quantile( subset( dfPlot, Group == x )$expA ); 
		return( ( qt[ "75%" ] - qt[ "25%" ] ) * 1.5 + qt[ "75%" ] ); } ) ) );
	
	yTick <- function(y) { 
		
		s <- 10^floor( log10( y ) ); 
		Exp_tick <- 0; 
		Exp_tick0 <- seq( 0, ceiling( y / s ) ) * s; 
		Exp_tick1 <- seq( 0, ceiling( y / s / 4 ) ) * s * 4; 
		Exp_tick2 <- seq( 0, ceiling( y / s / 5 ) ) * s * 5; 
		if( length( Exp_tick0 ) %in% c( 3, 4, 5 ) ) { Exp_tick <- Exp_tick0; return( Exp_tick ); }
		if( length( Exp_tick1 ) %in% c( 3, 4, 5 ) ) { Exp_tick <- Exp_tick1; return( Exp_tick ); }
		if( length( Exp_tick2 ) %in% c( 3, 4, 5 ) ) { Exp_tick <- Exp_tick2; return( Exp_tick ); }
		return( Exp_tick ); 
		
	}
	Exp_tick <- yTick( ymax ); 	
	
	if( length( Exp_tick ) == 1 ) { 
		
		s <- 10^( floor( log10( ymax ) ) - 1 ); 
		Exp_tick <- 0; 
		Exp_tick0 <- seq( 0, ceiling( ymax / s ) ) * s; 
		Exp_tick1 <- seq( 0, ceiling( ymax / s / 4 ) ) * s * 4; 
		Exp_tick2 <- seq( 0, ceiling( ymax / s / 5 ) ) * s * 5; 
		if( length( Exp_tick0 ) == 4 ) { Exp_tick <- Exp_tick0; }
		if( length( Exp_tick0 ) == 5 ) { Exp_tick <- Exp_tick0; }
		if( length( Exp_tick1 ) == 4 ) { Exp_tick <- Exp_tick1; }
		if( length( Exp_tick1 ) == 5 ) { Exp_tick <- Exp_tick1; }
		if( length( Exp_tick2 ) == 4 ) { Exp_tick <- Exp_tick2; }
		if( length( Exp_tick2 ) == 5 ) { Exp_tick <- Exp_tick2; }
		
	}
	print( ymax ); 
	print( Exp_tick ); 
	
	qt <- do.call( rbind, lapply( c( tlst0, tlst1 ), function(t) { 
		
		qt <- quantile( subset( dfPlot, Group == t )$expA ); 
		return( cbind( data.frame( Tissue=t ), data.frame( t( as.data.frame( qt ) ), check.names=FALSE ) ) ); 
		
		} ) )
	qt$Tissue <- unlist( lapply( as.character( qt$Tissue ), function(y) { return( gsub( "^NormalGTEX, ", "GTEx-", y ) ); } ) ); 
	write.table( qt, file=paste( fn, ".txt", sep="" ), append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 
	
}
plotEg0( "CLDN6" )
plotEg0( "CEACAM5" )
plotEg0( "CD276" )
	# --> Figure 3e

plotEg0( "MSLN" )
plotEg0( "TM4SF4" )
plotEg0( "GPC3" )
plotEg0( "TM4SF5" )
	# --> Figure 4c

plotEg0( "CA9" )
plotEg0( "SLC26A9" )
plotEg0( "GPA33" )
plotEg0( "TMIGD1" )
	# --> Figure 4g

plotEg0( "CD70" )
plotEg0( "CD27" )
plotEg0( "ULBP2" )
plotEg0( "KLRK1" )
	# --> Figure 6i







