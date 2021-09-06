#Figure 7a

library( "data.table" ); 
library( "reshape2" ); 
library( "ggplot2" ); 
library( "ggnewscale" ); 
library( "Hmisc" ); 
library( "ggrepel" ); 
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') { # Function to plot color bar
    scale = (length(lut)-1)/(max-min)

    dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
     y = (i-1)/scale + min
     rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}
library( "ggdendro" )
library( "gridExtra" ); 


# load LTMG stratified by cell type
flst <- system( "ls ./Summary |grep txt$", intern=TRUE )
t <- lapply( flst, function(fn) { 
	
	cht <- gsub( ".txt", "", fn ); 
	print( cht ); 
	
	scExp <- fread( file.path( "./Summary", fn ), data.table=FALSE, skip=0 ); 
	
	indC <- which( ! grepl( "Mean_", as.vector( t( scExp[ 1, ] ) ) ) & ! grepl( "ENSG", as.vector( t( scExp[ 2, ] ) ) ) ); 
	t <- scExp[ 1:2, indC ]; 
	cells <- data.frame( Num=as.numeric( as.character( as.vector( t( t[ 1, ] ) ) ) ) ); 
	row.names( cells ) <- as.vector( t( t[ 2, ] ) ); 
	
	t <- scExp[ 3:nrow( scExp ), indC ]; 
	
	ind <- which( cells[ , "Num" ] >= 50 ); # Tumors and non-malignant cell types containing <50 cells were excluded from the downstream analysis.
	r <- do.call( cbind, lapply( ind, function(x) { 
		
		r <- data.frame( t=as.numeric( as.character( t[ , x ] ) ) / cells[ x, "Num" ] ); 
		colnames( r ) <- row.names( cells )[ x ]; 
		return( r ); 
		
		} ) )
	row.names( r ) <- as.character( scExp[ 3:nrow( scExp ), 1 ] ); 
	print( r[ "ENSG00000153563", ] ); # CD8A
	
	Percdf <- melt( cbind( data.frame( GeneAcc=row.names( r ) ), r ), id="GeneAcc" ); 
	Percdf$cohort <- cht; 
	Percdf$variable <- unlist( lapply( as.character( Percdf$variable ), function(x) { return( ifelse( x == cht, "Cancer", x ) ); } ) ); 
	
	return( Percdf ); 
	
	} )
Percdf <- do.call( rbind, t ); 

cohorts <- unique( as.character( Percdf$cohort ) ); 


# calculate correlation based on cohort::cell::Frac

# define immune-related GESPs
dfAnno <- read.table( "immune_surface.V29.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
IA_GESP<- as.character( subset( dfAnno, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` != 0 )$`Ensembl ID` ); 
print( length( IA_GESP ) ) # 614

dfIA <- subset( Percdf, GeneAcc %in% IA_GESP );  
dfHIA <- acast( dfIA, GeneAcc ~ cohort + variable )
res2 <- rcorr( as.matrix( dfHIA ), type="spearman" ); 
dfRho <- melt( cbind( data.frame( Cell1=row.names( res2$r ) ), res2$r ), id="Cell1" ); 
colnames( dfRho ) <- c( "Cell1", "Cell2", "rho" ); 

t1 <- data.frame( Cohort1=unlist( lapply( as.character( dfRho$Cell1 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ), CellType1=unlist( lapply( as.character( dfRho$Cell1 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ) )
t2 <- data.frame( Cohort2=unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ), CellType2=unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ) )
dfRho <- cbind( dfRho, t1, t2 ); 
dfRho$CellTypeRank1 <- match( as.character( dfRho$CellType1 ), c( "Cancer", "Endothelial", "Fibroblast", "Adipocyte", "B cell", "T CD4", "T reg", "T CD8", "NK cell", "Granulocyte", "DC", "Macrophage", "Monocyte" ) )


plot1 <- function() { # plot heatmaps separately
	
	# remove HNSC and SKCM

	dfr <- res2$r; 
	print( head( row.names( dfr ) ) ); 
	ind <- which( ! unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ) %in% c( "HNSC", "SKCM" ) ); 
	dfr <- dfr[ ind, ind ]; 
	
	dfCohortCell <- data.frame( CohortCell=row.names( dfr ), 
		Cohort=unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ), 
		Cell=unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ) )
	
	tRho <- data.frame(); 
	t <- lapply( c( "Cancer", "Endothelial", "Fibroblast", "Adipocyte", "B cell", "T CD4", "T reg", "T CD8", "NK cell", "Granulocyte", "DC", "Macrophage", "Monocyte", 
		"Macrophage/Monocyte", "Macrophage/Monocyte/DC", "Monocyte/DC" ), function(x) { 
		
		if( x == "Macrophage/Monocyte" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Macrophage", "Monocyte" ) ); 
		} else if( x == "Macrophage/Monocyte/DC" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Macrophage", "Monocyte", "DC" ) ); 
		} else if( x == "Monocyte/DC" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Monocyte", "DC" ) ); 
		} else if( x == "T cells" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "T CD4", "T reg", "T CD8" ) ); 
		} else { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% x ); 
		}
		
		dfr <- dfr[ ind, ind ]; 
		print( round( quantile( as.vector( as.matrix( dfr ) ) ), 2 ) ); 
		
		hc <- hclust( dist( t( dfr ) ), method="complete" ); 
		clo <- colnames( dfr )[ hc$order ]; 
		print( clo ); 
		
		dhc <- as.dendrogram( hc )
		ddata <- dendro_data( dhc, type="rectangle" )
		dfSeg <- segment( ddata ); 
		
		dfRho <- melt( cbind( data.frame( Cell1=row.names( dfr ) ), dfr ), id="Cell1" ); 
		colnames( dfRho ) <- c( "Cell1", "Cell2", "rho" ); 
		dfRho$CellRank1 <- match( as.character( dfRho$Cell1 ), clo )
		dfRho$CellRank2 <- match( as.character( dfRho$Cell2 ), clo )
		dfRho$CellType1 <- unlist( lapply( as.character( dfRho$Cell1 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ); 
		dfRho$CellType2 <- unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ); 
		dfRho$CellCohort1 <- unlist( lapply( as.character( dfRho$Cell1 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ); 
		dfRho$CellCohort2 <- unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ); 
		
		p1 <- ggplot() + 
			geom_rect( data=dfRho, mapping=aes( xmin=CellRank1 - 0.5, xmax=CellRank1 + 0.5, ymin=CellRank2 - 0.5, ymax=CellRank2 + 0.5, fill=as.numeric( as.character( rho ) ) ), color=NA, alpha=1 ) +
			scale_fill_gradientn( guide=guide_colourbar(), name="Rho", colours=
				c( colorRampPalette( c( "#0000FF", "#1E90FF" ) )( 60 ), colorRampPalette( c( "#1E90FF", "orange" ) )( 10 ), colorRampPalette( c( "orange", "#FF3333" ) )( 10 ), colorRampPalette( c( "#FF3333", "red" ) )( 20 ) ), 
				na.value="#E5E5E5", limits=c( 0, 1 ) ) +
			theme( legend.position="none", legend.justification="left", aspect.ratio=1 / 1.6,
				plot.title=element_text( size=5 ), 
				panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(),
				axis.line.x=element_blank(), axis.line.y=element_blank(), 
				axis.text.x=element_blank(), axis.text.y=element_blank(), 
				axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
				plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
			xlab( "" ) + ylab( "" ) + labs( title=x )
		# p1; dev.off(); 
		
		if( length( unique( c( as.character( dfRho$CellType1 ), as.character( dfRho$CellType2 ) ) ) ) == 1 ) { # single immune cell type
			
			out <- dfRho; 
			out$`Cell type` <- x; 
			
			t <- sort( unique( c( as.character( dfRho$CellCohort1 ), as.character( dfRho$CellCohort2 ) ) ) ); 
			out$t1 <- match( as.character( dfRho$CellCohort1 ), t ); 
			out$t2 <- match( as.character( dfRho$CellCohort2 ), t ); 
			out <- subset( out, t1 < t2 ); 
			
			out$`Cancer 1` <- as.character( out$CellCohort1 ); 
			out$`Cancer 2` <- as.character( out$CellCohort2 ); 
			out$`Correlation coefficient` <- as.character( out$rho ); 
			
			# print( showd( out ) ); 
			write.table( out[ , c( "Cell type", "Cancer 1", "Cancer 2", "Correlation coefficient" ) ], 
				file=paste( "t", gsub( "/", "_", x ), "Rho.txt" ), append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, na="", sep="\t" ); 
			
		} else { # mixed immune cell types
			
			out <- dfRho; 
			
			t <- sort( unique( c( as.character( dfRho$CellCohort1 ), as.character( dfRho$CellCohort2 ) ) ) ); 
			out$t1 <- match( as.character( dfRho$CellCohort1 ), t ); 
			out$t2 <- match( as.character( dfRho$CellCohort2 ), t ); 
			out <- subset( out, t1 < t2 ); 
			
			out$`Cancer 1` <- as.character( out$CellCohort1 ); 
			out$`Cancer 2` <- as.character( out$CellCohort2 ); 
			out$`Correlation coefficient` <- as.character( out$rho ); 
			
			out$`Cell 1` <- as.character( out$CellType1 ); 
			out$`Cell 2` <- as.character( out$CellType2 ); 
			
			# print( showd( out ) ); 
			write.table( out[ , c( "Cancer 1", "Cell 1", "Cancer 2", "Cell 2", "Correlation coefficient" ) ], 
				file=paste( "t", gsub( "/", "_", x ), "Rho.txt" ), append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, na="", sep="\t" ); 
			
		}
		
		tRho <<- rbind( tRho, cbind( data.frame( Test=x ), dfRho ) ); 
		return( p1 ); 
		
		} )
	
	pdf( "tmp.pdf", width=12, height=9.6 ); 
	do.call( grid.arrange, t ) 
	dev.off()
	
	write.table( tRho, file="tmp.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 
		# -> Figure 7A
	
}
plot1()

plot1 <- function() { # plot clustering tree separately
	
	# remove HNSC and SKCM

	dfr <- res2$r; 
	print( head( row.names( dfr ) ) ); 
	ind <- which( ! unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ) %in% c( "HNSC", "SKCM" ) ); 
	dfr <- dfr[ ind, ind ]; 
	
	dfCohortCell <- data.frame( CohortCell=row.names( dfr ), 
		Cohort=unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ), 
		Cell=unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ) )
	
	t <- lapply( c( "Cancer", "Endothelial", "Fibroblast", "Adipocyte", "B cell", "T CD4", "T reg", "T CD8", "NK cell", "Granulocyte", "DC", "Macrophage", "Monocyte", 
		"Macrophage/Monocyte", "Macrophage/Monocyte/DC", "Monocyte/DC" ), function(x) { 
		
		if( x == "Macrophage/Monocyte" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Macrophage", "Monocyte" ) ); 
		} else if( x == "Macrophage/Monocyte/DC" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Macrophage", "Monocyte", "DC" ) ); 
		} else if( x == "Monocyte/DC" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Monocyte", "DC" ) ); 
		} else if( x == "T cells" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "T CD4", "T reg", "T CD8" ) ); 
		} else { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% x ); 
		}
		
		dfr <- dfr[ ind, ind ]; 
		print( round( quantile( as.vector( as.matrix( dfr ) ) ), 2 ) ); 
		
		hc <- hclust( dist( t( dfr ) ), method="complete" ); 
		clo <- colnames( dfr )[ hc$order ]; 
		print( clo ); 
		print( rev( clo ) ); 
		
		dhc <- as.dendrogram( hc )
		ddata <- dendro_data( dhc, type="rectangle" )
		dfSeg <- segment( ddata ); 
		
		dfRho <- melt( cbind( data.frame( Cell1=row.names( dfr ) ), dfr ), id="Cell1" ); 
		colnames( dfRho ) <- c( "Cell1", "Cell2", "rho" ); 
		dfRho$CellRank1 <- match( as.character( dfRho$Cell1 ), clo )
		dfRho$CellRank2 <- match( as.character( dfRho$Cell2 ), clo )
		dfRho$CellType1 <- unlist( lapply( as.character( dfRho$Cell1 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ); 
		dfRho$CellType2 <- unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ); 
		dfRho$CellCohort2 <- unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ); 
		
		p1 <- ggplot() + 
			geom_segment( data=dfSeg, aes( x=y, y=x, xend=yend, yend=xend ) ) + 
			theme( legend.position="none", legend.justification="left", aspect.ratio=34 / 4.8,
				panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(),
				axis.line.x=element_blank(), axis.line.y=element_blank(), 
				axis.text.x=element_blank(), axis.text.y=element_blank(), 
				axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
				plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
			xlab( "" ) + ylab( "" ) + labs( title="" )
		# p1; dev.off(); 
		
		return( p1 ); 
		
		} )
	
	pdf( "tmp.pdf", width=10, height=9.6 ); 
	do.call( grid.arrange, t ) 
	dev.off()
	
	pdf( "tmp.pdf", width=10, height=20 ); 
	do.call( grid.arrange, t ) 
	dev.off()
	
}
plot1()

plot1 <- function() { # plot violin
	
	# remove HNSC and SKCM

	dfr <- res2$r; 
	print( head( row.names( dfr ) ) ); 
	ind <- which( ! unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ) %in% c( "HNSC", "SKCM" ) ); 
	dfr <- dfr[ ind, ind ]; 
	
	dfCohortCell <- data.frame( CohortCell=row.names( dfr ), 
		Cohort=unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ), 
		Cell=unlist( lapply( row.names( dfr ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ) )
	print( length( sort( unique( as.character( dfCohortCell$Cohort ) ) ) ) ); # 13
	
	t <- lapply( c( "Cancer", "Endothelial", "Fibroblast", "Adipocyte", "B cell", "T CD4", "T reg", "T CD8", "NK cell", "Granulocyte", "DC", "Macrophage", "Monocyte" ), function(x) { 
		
		if( x == "Macrophage/Monocyte" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Macrophage", "Monocyte" ) ); 
		} else if( x == "Macrophage/Monocyte/DC" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Macrophage", "Monocyte", "DC" ) ); 
		} else if( x == "Monocyte/DC" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "Monocyte", "DC" ) ); 
		} else if( x == "T cells" ) { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% c( "T CD4", "T reg", "T CD8" ) ); 
		} else { 
			ind <- which( as.character( dfCohortCell$Cell ) %in% x ); 
		}
		
		dfr <- dfr[ ind, ind ]; 
		print( round( quantile( as.vector( as.matrix( dfr ) ) ), 2 ) ); 
		
		dfRho <- melt( cbind( data.frame( Cell1=row.names( dfr ) ), dfr ), id="Cell1" ); 
		colnames( dfRho ) <- c( "Cell1", "Cell2", "rho" ); 

		t1 <- data.frame( Cohort1=unlist( lapply( as.character( dfRho$Cell1 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ), CellType1=unlist( lapply( as.character( dfRho$Cell1 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ) )
		t2 <- data.frame( Cohort2=unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 1 ] ); } ) ), CellType2=unlist( lapply( as.character( dfRho$Cell2 ), function(x) { return( unlist( strsplit( x, "_" ) )[ 2 ] ); } ) ) )
		dfRho <- cbind( dfRho, t1, t2 ); 
		dfRho$CohortRank1 <- match( as.character( dfRho$Cohort1 ), sort( unique( c( as.character( dfRho$Cohort1 ), as.character( dfRho$Cohort2 ) ) ) ) )
		dfRho$CohortRank2 <- match( as.character( dfRho$Cohort2 ), sort( unique( c( as.character( dfRho$Cohort1 ), as.character( dfRho$Cohort2 ) ) ) ) )
		
		return( dfRho ); 
		
		} )
	dfRho <- do.call( rbind, t ); 
	
	t <- c( "Cancer", "Adipocyte", "Endothelial", "Fibroblast", "Granulocyte", "NK cell", "B cell", "T CD8", "T CD4", "T reg", "Monocyte", "Macrophage", "DC" ); 
	names( t ) <- seq( 1, length( t ) ); 
	dfRho$CellType1Rank <- match( as.character( dfRho$CellType1 ), t ); 
	
	write.table( data.frame( x="" ), file="tmp.txt", append=FALSE, quote=F, row.names=FALSE, col.names=FALSE, na="", sep="\t" ); 
		# --> Figure 7D
	x <- lapply( t, function(x) { 
		
		write.table( data.frame( t( rbind( data.frame( x=x ), data.frame( x=as.character( subset( dfRho, CellType1 == x & CellType2 == x & CohortRank1 < CohortRank2 )$rho ) ) ) ), check.names=FALSE ), 
			file="tmp.txt", append=TRUE, quote=F, row.names=FALSE, col.names=FALSE, na="", sep="\t" ); 
		
		} ) 
	
	g <- ggplot( data=subset( dfRho, CellType1 == CellType2 & CohortRank1 < CohortRank2 ), aes( x=CellType1Rank, y=rho, fill=CellType1 ) ) +
		geom_violin( trim=F, size=0.2, show.legend =F, width=1.1, colour="#FFFFFF", alpha=0.84 ) + 
		geom_jitter( width=0.25, colour="#BFBFBF", size=0.4 ) + 
		stat_summary( fun.y=median, geom="point", size=0.6, color="black")+
		scale_fill_manual( name="Cell type", values=c( "Cancer"="#FF0000", "Adipocyte"="#cf6e5d", "B cell"="#cd6785", "DC"="#563687", "Endothelial"="#ab5fa6", "Fibroblast"="#C15064", "Granulocyte"="#86271a", "Macrophage"="#d49830", "Monocyte"="#4cc0a9", "NK cell"="#666eb3", "T CD4"="#85ab3f", "T CD8"="#c0a953", "T reg"="#7dbc5e" ) )	+
		scale_x_continuous( breaks=seq( 1, length( t ) ), labels=NULL, limits=c( 0.4, length( t ) + 0.6 ), expand=c( 0, 0 ) ) +
		labs( y="", x="", title="" ) + 
		theme_classic() + 
		theme( legend.position="none", aspect.ratio=48/287, 
			axis.line.x=element_line( colour="black", size=0.8, linetype=1 ), axis.line.y=element_line( colour="black", size=0.8, linetype=1 ), 
			panel.border=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
			axis.ticks.x=element_line( colour="black", size=0.8, linetype=1 ), axis.ticks.y=element_line( colour="black", size=0.8, linetype=1 ), axis.ticks.length=unit( 0.2, "cm" ), 
			axis.text.x=element_text( colour="black", size=6, angle=45, hjust=1, vjust=0.5 ), axis.text.y=element_text( colour="black", size=6 ) )
	ggsave( file="tmp.svg", plot=g, width=9, height=5 )
	
}
plot1()

#Figure 7e
library( "ggplot2" )


# define immune-related GESPs
dfAnno <- read.table( "immune_surface.V29.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
IA_GESP<- as.character( subset( dfAnno, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` != 0 )$`Ensembl ID` ); 
print( length( IA_GESP ) )
# 614


PC1 <- read.table( "PC1.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
# 58676 genes 

PC1$GeneAcc <- unlist( lapply( as.character( PC1$Gene ), function(x) { xx <- unlist( strsplit( x, " \\(" ) )[ 2 ]; return( unlist( strsplit( xx, "\\)" ) )[ 1 ] ); } ) )
PC1 <- cbind( PC1, dfAnno[ match( as.character( PC1$GeneAcc ), as.character( dfAnno$`Ensembl ID` ) ), c( "Immunological Accessory Molecules", "Surfaceome" ) ] )


tab <- table( subset( PC1, GeneAcc %in% IA_GESP )$GroupNew ); 
print( tab ); 
#     Bimodal-Like Expressed-In-All Lineage-Enriched    Not-Expressed            Other     Right-Skewed 
#               55               44              116              126              223               50
print( round( tab / sum( tab ) * 100, 1 ) ); 
#     Bimodal-Like Expressed-In-All Lineage-Enriched    Not-Expressed            Other     Right-Skewed 
#              9.0              7.2             18.9             20.5             36.3              8.1 

print( as.data.frame( tab ) ); 
print( as.data.frame( round( tab / sum( tab ), 3 ) ) ); 
print( sum( tab[ c( "Expressed-In-All", "Bimodal-Like", "Lineage-Enriched", "Other", "Right-Skewed" ) ] ) ); # [1] 488
print( sum( tab[ c( "Expressed-In-All", "Bimodal-Like", "Lineage-Enriched", "Other", "Right-Skewed" ) ] ) / sum( tab ) );  # [1] 0.7947883

print( sum( tab[ c( "Bimodal-Like", "Lineage-Enriched", "Other", "Right-Skewed" ) ] ) ); # [1] 444
print( sum( tab[ c( "Bimodal-Like", "Lineage-Enriched", "Other", "Right-Skewed" ) ] ) / sum( tab ) );  # [1] 0.723127

f <- tab[ "Not-Expressed" ] / sum( tab ); 
dfP <- data.frame( TEgroup="Not-Expressed", xmin=0, xmax=f, ymin=0, ymax=1 ); 

tab <- table( as.character( subset( PC1, GeneAcc %in% IA_GESP & GroupNew != "Not-Expressed" )$GroupNew ) ); 
go <- rev( c( "Expressed-In-All", "Lineage-Enriched", "Right-Skewed", "Bimodal-Like", "Other" ) ); 
Pct <- as.data.frame( tab[ go ] / sum( tab[ go ] ) )
colnames( Pct ) <- c( "GeneClass", "value" ); 
Pct <- do.call( rbind, lapply( seq( 1, nrow( Pct ) ), function(y) { 
	
	ymin <- ifelse( y == 1, 0, sum( Pct[ 1:(y-1), ]$value ) ); 
	ymax <- ymin + Pct[ y, ]$value; 
	return( cbind( Pct[ y, ], data.frame( ymin=ymin, ymax=ymax ) ) ); 
	
	} ) ); 
dfP <- rbind( dfP, data.frame( TEgroup=as.character( Pct$GeneClass ), xmin=f, xmax=1, ymin=Pct$ymin, ymax=Pct$ymax ) ); 

p1 <- ggplot( ) +
	geom_rect( data=dfP, aes( fill=as.character( TEgroup ), xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax ), color="#FFFFFF" ) +
	scale_x_continuous( limits=c( 0, 1 ), breaks=NULL ) +
	scale_y_continuous( limits=c( 0, 1 ), breaks=NULL, labels=NULL ) +
	xlab( "" ) + ylab( "" ) + ggtitle( paste( "Immune GESPs (n=", sum( tab ), ")", sep="" ) ) +
	scale_fill_manual( name="", values=c( "Expressed-In-All"="#f08585", "Bimodal-Like"="#f6a104", "Right-Skewed"="#f8b437", "Other"="#9c6d17", "Lineage-Enriched"="#facd79", "Not-Expressed"="#b3b1b1" ) ) + 
	guides( fill=guide_legend( nrow=2, byrow=TRUE ) ) + 
	theme( legend.position="bottom", aspect.ratio=1, 
		panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
		panel.border=element_blank(), panel.background=element_blank(), 
		axis.line.x=element_line( size=0.8, colour="black" ), axis.line.y=element_line( size=0.8, colour="black" ), 
		axis.text.x=element_text( size=8 ), axis.text.y=element_text( size=8 ), 
		axis.ticks.x=element_line( colour="black", size=0.5, linetype=1 ), axis.ticks.y=element_line( colour="black", size=0.5, linetype=1 ), axis.ticks.length=unit( 0.12, "cm" ), 
		plot.margin=unit( c( 0, 0, 0, 1 ), "cm" ) )
ggsave( file="tmp1.svg", plot=p1, width=6, height=7 )
# 	--> Figure 7e

t <- as.character( subset( PC1, GeneAcc %in% IA_GESP )$GroupNew ); 
t <- factor( t, levels=c( "Expressed-In-All", "Not-Expressed", "Lineage-Enriched", "Right-Skewed", "Bimodal-Like", "Other" ) ); 
write.table( cbind( as.data.frame( table( t ) ), as.data.frame( table( t ) / length( t ) ) ), file="tmp2.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 
	# -> Figure 7e



