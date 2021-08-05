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

#Figure 7g
library( "miceadds" );
library( "DescTools" )
library( "limma" )
library( "reshape2" );
library( "ggplot2" ); 
library( "ggnewscale" ); 


sample_info <- read.csv( "sample_info.csv", header=TRUE, check.names=FALSE ); 

geneinfo <- read.table( "geneinfo.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 

Hallmark <- read.table( "CellSyst.2015.MSigDB.Hallmark.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 


# load ISG gene list 
ISG.gene <- as.vector( as.matrix( read.table( "ISG.gene.lst", header=FALSE, check.names=FALSE, sep="\t", quote="", comment.char="" ) ) );

# load HALLMARK gene sets
msigdb <- read.table( "msigdb.v6.0.symbols.gmt", sep="\t", fill=TRUE, header=F, 
	col.names=c( "GeneSet", "Link", 1 : max(count.fields( "msigdb.v6.0.symbols.gmt", sep = "\t") ) ) )
msigdb <- subset( msigdb, grepl( "^HALLMARK_", GeneSet ) ); 
row.names( msigdb ) <- as.character( msigdb$GeneSet ); 

t <- apply( msigdb[ , 3 : ncol( msigdb ) ], 1, function(x) { 
	
	x <- unique( sort( as.character( x ) ) ); 
	return( setdiff( x, "" ) ); 
	
	} )
names( t ) <- row.names( msigdb ); 

t$ISG <- ISG.gene; 
geneSets <- t; 


# load expression from DepMap
load.Rdata( filename="DepMap_ExpFull.RData", "Exp" );
dfS <- cbind( Exp$dfS, data.frame( Tissue=unlist( lapply( as.character( Exp$dfS$CCLE_Name ), function(x) { 
	
	xx <- unlist( strsplit( x, "_" ) ); 
	return( paste( xx[ 2 : length( xx ) ], collapse="_" ) ); 
	
	} ) ) ) )

# manually correct sample list included for downstream analysis 
sample <- read.table( "sample.Exp.LZ.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
dfS$TissueCorrected <- as.character( sample[ match( as.character( dfS$DepMap_ID ), as.character( sample$DepMap_ID ) ), "L Zhang Cancer Types" ] )

dfS$Tissue <- as.character( dfS$TissueCorrected ); 
tab <- table( sort( as.character( dfS$Tissue ) ) ); 

# select cell lines to include 
ind <- which( as.character( dfS$TissueCorrected ) %in% setdiff( names( tab[ tab >= 3 ] ), c( "FIBROBLAST", "" ) ) ); 
df <- Exp$df[ ind, ]; 
dfS <- dfS[ ind, ]; 
tab <- table( sort( as.character( dfS$Tissue ) ) ); 

geneRow <- Exp$dfG; 


# calculate signatures of gene sets
t <- lapply( geneSets, function(gs) { 
	
	indR <- match( gs, as.character( geneRow$GeneName ) ); 
	indR <- indR[ ! is.na( indR ) ]; 
	Exp <- df[ , indR ]; 
	geneRow <- geneRow[ indR, ]; 
	colnames( Exp ) <- as.character( geneRow$GeneName ); 

	geneExpZ <- apply( Exp, 2, function(x) { 
		
		mad <- MeanAD( x ); 
		med <- median( x ); 
		return( ( x - med ) / ( 1.253314 * mad ) ); 
		
		} )
	t <- data.frame( t=unlist( apply( geneExpZ, 1, mean ) ) ); 
	return( t ); 
	
	} )
t <- do.call( cbind, t ); 
row.names( t ) <- row.names( df ); 
colnames( t ) <- names( geneSets ); 

sigs <- t; 


# association test 
t <- lapply( colnames( sigs ), function( gs ) { 
	
	sml <- sigs[ , gs ]; 

	ind <- which( ! is.na( sml ) ); 
	sml <- sml[ ind ]; 

	design <- model.matrix( ~ sml ); 
	
	fit <- lmFit( data.frame( t( df[ ind, ] ) ), design ); 
	fit <- eBayes( fit, robust=TRUE ); 
	tT <- subset( topTable( fit, adjust="BH", sort.by="p", number=ncol( df ) ), select=c( "adj.P.Val","P.Value","t","B","logFC" ) ); 
	tT <- cbind( data.frame( GeneSet=gs, Gene=row.names( tT ) ), geneRow[ match( row.names( tT ), colnames( df ) ), ], tT ); 
	
	return( tT ); 
	
	} )
tT <- do.call( rbind, t ); 

t <- apply( tT, 1, function(x) { 
	
	if( as.character( x[ "GeneName" ] ) %in% geneSets[[ as.character( x[ "GeneSet" ] ) ]] ) { return( 1 ); 
	} else { return( 0 ); }
	
	} )
tT$Include <- unlist( t ); 

tC <- cor( df, sigs )

t1 <- apply( tT, 1, function(x) { return( paste( as.character( x[ "Gene" ] ), as.character( x[ "GeneSet" ] ), sep="::" ) ); } )
t2df <- melt( cbind( data.frame( Gene=row.names( tC ) ), tC ), id="Gene" ); 
t2 <- apply( t2df, 1, function(x) { return( paste( as.character( x[ "Gene" ] ), as.character( x[ "variable" ] ), sep="::" ) ); } )
tT$r <- t2df[ match( t1, t2 ), "value" ]; 


# clustering using only expressed IAMs 

TE6 <- read.table( "CCLE.TE.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 

tmp <- as.character( subset( TE6, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` == 1 & Group %in% c( "Expressed-In-All", "Bimodal-Like", "Lineage-Enriched", "Other", "Right-Skewed" ) )$GeneAcc )
dfP <- subset( tT, EnsemblID %in% tmp ); 
dfPr <- acast( dfP, GeneName ~ GeneSet, value.var="r" )


hc <- hclust( dist( t( dfPr ) ), method="complete" ); 
ho <- colnames( dfPr )[ hc$order ]; 
hg <- cutree( hc, k=7 )


# adjust groups / orders of hallmarks 

# ISG as a last outlier 
hg[ "ISG" ] <- 8; 
ho <- c( setdiff( ho, "ISG" ), "ISG" ); 
 
# merge ccluster#2 and cluster#3
hg[ names( which( hg == 3 ) ) ] <- 2; 

# shift two (2) pathways from cluster#1 to cluster#5
hg[ ho[ head( sort( match( names( which( hg == 1 ) ), ho ) ), n=2 ) ] ] <- 5; 

# re-number clusters 
print( table( hg ) ); 
t <- names( hg ); 
hg <- c( 6, 2, 4, 5, 3, 1, 7 )[ match( hg, c( 1, 2, 4, 5, 6, 7, 8 ) ) ];
names( hg ) <- t; 
print( table( hg ) ); 


hc <- hclust( dist( as.matrix( dfPr ) ), method="complete" ); 
go <- row.names( dfPr )[ hc$order ]; 
gg <- cutree( hc, k=8 )


# adjust groups / orders of genes 

# merge ccluster#6 and cluster#3
gg[ names( which( gg == 6 ) ) ] <- 3; 

# merge ccluster#1 and cluster#2
gg[ names( which( gg == 2 ) ) ] <- 1; 

# merge ccluster#4 and cluster#8
gg[ names( which( gg == 8 ) ) ] <- 4; 

# re-number clusters 
print( table( gg ) ); 
t <- names( gg ); 
gg <- c( 2, 1, 3, 5, 4 )[ match( gg, c( 1, 3, 4, 5, 7 ) ) ];
names( gg ) <- t; 
print( table( gg ) ); 

# reverse gene order
go <- rev( go ); 
gg <- rev( gg ); 


print( dim( dfPr ) ); 
print( head( dfPr[ , 1:4 ] ) )

print( head( go ) ); 
print( head( gg ) ); 
print( head( gg[ go ] ) ); 
print( head( dfPr[ rev( go ), 1:4 ] ) )

print( head( ho ) ); 
print( head( hg ) ); 
print( head( hg[ ho ] ) ); 

write.table( dfPr[ rev( go ), ho ], file="tmp1.txt", append=FALSE, quote=F, row.names=TRUE, col.names=NA, na="", sep="\t" ); 
# 	--> Supplementary Data 31. Correlation between expression levels of membrane-bound IAMs and 50 “hallmark” gene sets in cancer cell lines

MasterTableV31 <- read.table( "MasterTableV31.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
t <- rev( go ); 
write.table( MasterTableV31[ match( t, MasterTableV31$`Gene Name` ), ], file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 

general <- read.table( "general.txt", header=FALSE, check.names=FALSE, sep="\t", quote="", comment.char="" ); 
t <- rev( go ); 
write.table( general[ match( t, general$`V2` ), ], file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 

write.table( data.frame( Gene=rev( go ), Cluster=c( "A", "B", "C", "D", "E" )[ match( gg[ rev( go ) ], c( 1, 2, 3, 4, 5 ) ) ] ), file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 

write.table( data.frame( Pathway=ho, Cluster=hg[ ho ], 
	Category=as.character( Hallmark[ match( gsub( "^HALLMARK_", "", ho ), as.character( Hallmark$`Hallmark Name` ) ), "Process Category" ] ) ), 
	file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 


dfAnnoImmune <- read.table( "immune_surface.ImmuneV29.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 


dfRho <- melt( cbind( data.frame( Gene=row.names( dfPr ) ), dfPr ), id="Gene" ); 
colnames( dfRho ) <- c( "Gene", "HALL", "r" ); 
dfRho$GeneRank <- match( as.character( dfRho$Gene ), go )
dfRho$HALLRank <- match( as.character( dfRho$HALL ), ho )


marginRow <- nrow( dfPr ) / 60 / ( 115 / 356 ); 
marginCol <- ncol( dfPr ) / 60; 
annoStartRow <- ( nrow( dfPr ) / 60 / ( 115 / 356 ) ) * (-1);
annoStartCol <- ( ncol( dfPr ) / 60 ) * (-1);


p1 <- ggplot() + 
	geom_rect( data=dfRho, mapping=aes( xmin=HALLRank - 0.5, xmax=HALLRank + 0.5, ymin=GeneRank - 0.5, ymax=GeneRank + 0.5, fill=as.numeric( as.character( r ) ) ), color=NA, alpha=1 ) +
	scale_fill_gradientn( guide=guide_colourbar( barwidth=0.5, barheight=3 ), name="Rho", colours=
		c( colorRampPalette( c( "blue", "#1E90FF" ) )( 35 ), colorRampPalette( c( "#8080FF", "#FFFFFF" ) )( 15 ), colorRampPalette( c( "#FFFFFF", "#FF8080" ) )( 15 ), colorRampPalette( c( "#FF8080", "#FF0000" ) )( 35 ) ), 
		na.value="#E5E5E5", limits=c( -1, 1 ) ) +

	new_scale_fill() +


	geom_rect( data=data.frame( HALL=ho, rank=seq( 1, length( ho ) ), group=paste( "h", hg[ ho ], sep="" ) ), 
		mapping=aes( xmin=rank - 0.5, xmax=rank + 0.5, ymin=annoStartRow, ymax=annoStartRow - ( marginRow * 1 ), fill=as.character( group ) ), color=NA, alpha=1 ) +
	geom_text( data=data.frame( y=annoStartRow - ( marginRow * 2 ), x=unlist( lapply( sort( unique( hg ) ), function(x) { return( median( match( names( which( hg == x ) ), ho ) ) ); } ) ), label=sort( unique( hg ) ) ), 
		aes( x, y, label=label ), hjust=0.5, vjust=0.5 ) +
		

	geom_rect( data=data.frame( Gene=go, rank=seq( 1, length( go ) ), group=paste( "g", gg[ go ], sep="" ) ), 
		mapping=aes( ymin=rank - 0.5, ymax=rank + 0.5, xmin=annoStartCol - ( marginCol * 0 ), xmax=annoStartCol - ( marginCol * 1 ), fill=as.character( group ) ), color=NA, alpha=1 ) +
	geom_text( data=data.frame( x=annoStartCol - ( marginCol * 2 ), y=unlist( lapply( sort( unique( gg ) ), function(x) { return( median( match( names( which( gg == x ) ), go ) ) ); } ) ), label=sort( unique( gg ) ) ), 
		aes( x, y, label=label ), hjust=0, vjust=0.5 ) +

	
	scale_fill_manual( guide=FALSE, name="Cluster", values=c( "h1"="#fb0015", "h2"="#CC008A", "h3"="#9d00ff", 
		"h4"="#0d00ff", "h5"="#0d6eff", "h6"="#22ffc4", "h7"="#22ff28", "g1"="#f6e41b", "g2"="#63ca45", "g3"="#209673", 
		"g4"="#245d7b", "g5"="#342770" ) ) + 

# 	geom_text( data=data.frame( x=seq( 1, length( ho ) ), y=nrow( dfPr ) + ( marginRow ), label=as.character( Hallmark[ match( gsub( "^HALLMARK_", "", ho ), as.character( Hallmark$`Hallmark Name` ) ), "Process Category" ] ) ), 
# 		aes( x, y, label=label ), hjust=0, vjust=0.5, angle=90, size=2 ) +

# 	scale_x_continuous( breaks=seq( 1, length( ho ) ), labels=gsub( "^HALLMARK_", "", ho ) )+
	scale_x_continuous( breaks=NULL )+
	scale_y_continuous( breaks=NULL, limits=c( annoStartRow - ( marginRow * 6 ), nrow( dfPr ) + 1 + ( nrow( dfPr ) / 10 )  ) )+
	theme( legend.position="right", legend.justification="left", 
		aspect.ratio=115 / 356,
		plot.title=element_text( size=12 ), 
		panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(),
		axis.line.x=element_blank(), axis.line.y=element_blank(), 
		axis.text.x=element_text( colour="black", size=4.8, angle=90, hjust=1, vjust=0.5 ), axis.text.y=element_blank(), 
		axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
		plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
	xlab( "HALLMARK" ) + ylab( paste0( "Gene (n=", nrow( dfPr ), ")" ) ) + labs( title=paste0( "CCLE (n=", nrow( df ), ")" ) )
ggsave( file="tmp1.svg", plot=p1, width=10, height=5 )

write.table( data.frame( Pathway=ho, Category=as.character( Hallmark[ match( gsub( "^HALLMARK_", "", ho ), as.character( Hallmark$`Hallmark Name` ) ), "Process Category" ] ) ), 
	file="tmp1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, na="", sep="\t" ); 


# calculate enrichment
# GeneClass <- c( "Adhesion Molecules", "Fc Receptors", "Tetraspanins", "Cytokine", "Cytokine Receptors", "Chemokine", "Chemokine Receptors", "Co-stimulator/Co-inhibitor", "Membrane-bound pattern recognition receptors" ); 
GeneClass <- c( "Adhesion Molecules", "Fc Receptors", "Tetraspanins", "Cytokine/Cytokine Receptors", "Chemokine/Chemokine Receptors", "Co-stimulator/Co-inhibitor", "Membrane-bound pattern recognition receptors" ); 
ggo <- unique( gg[ go ] ); 
AllGeneNames <- as.character( subset( TE6, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` == 1 & Group %in% c( "Expressed-In-All", "Bimodal-Like", "Lineage-Enriched", "Other", "Right-Skewed" ) )$`Gene Name` )
t <- lapply( seq( 1, length( GeneClass ) ), function(i) { 
	
	geneClass <- GeneClass[ i ]; 
	
	if( geneClass == "Cytokine/Cytokine Receptors" ) { 
		genes <- unique( dfAnnoImmune[ as.character( dfAnnoImmune[ , "Cytokine" ] ) != "", "Gene Name" ], dfAnnoImmune[ as.character( dfAnnoImmune[ , "Cytokine Receptors" ] ) != "", "Gene Name" ] ); 
	} else if( geneClass == "Chemokine/Chemokine Receptors" ) { 
		genes <- unique( dfAnnoImmune[ as.character( dfAnnoImmune[ , "Chemokine" ] ) != "", "Gene Name" ], dfAnnoImmune[ as.character( dfAnnoImmune[ , "Chemokine Receptors" ] ) != "", "Gene Name" ] ); 
	} else { 
		genes <- dfAnnoImmune[ as.character( dfAnnoImmune[ , geneClass ] ) != "", "Gene Name" ]; 
	}
	
	return( do.call( rbind, lapply( ggo, function(cluster) { 
		
		tab <- table( AllGeneNames %in% names( which( gg == cluster ) ), AllGeneNames %in% genes ); 
		fit <- fisher.test( tab ); 
		return( data.frame( GeneClass=geneClass, Cluster=cluster, OR=as.numeric( as.character( fit$estimate ) ), pvalue=as.numeric( as.character( fit$p.value ) ) ) ); 
		
		} ) ) ); 
	
	} )
fit <- do.call( rbind, t ); 
fit$BH <- p.adjust( fit$pvalue, method="BH" ); 
fit$xRank <- match( fit$GeneClass, GeneClass ); 
fit$yRank <- match( fit$Cluster, ggo ); 
print( quantile( subset( fit, OR > 1 )$OR ) ); 
print( quantile( -log10( subset( fit, OR > 1 )$pvalue ) ) ); 
p4 <- ggplot() + 
	geom_point( data=subset( fit, OR > 1 ), aes( x=xRank, y=unlist( lapply( Cluster, function(x) { return( median( match( names( which( gg == x ) ), go ) ) ); } ) ), 
		size=OR, fill=-log10( pvalue ) ), color="#D9D9D9", shape=21 ) +
	geom_point( data=subset( fit, OR > 1 & BH < 0.1 ), aes( x=xRank, y=unlist( lapply( Cluster, function(x) { return( median( match( names( which( gg == x ) ), go ) ) ); } ) ), 
		size=OR, fill=-log10( pvalue ) ), shape=21, color="#000000" ) +
	scale_fill_gradientn( guide=guide_colourbar( barwidth=0.5, barheight=3 ), name="Significance", colours=
		c( colorRampPalette( c( "#FFFFFF", "#FF8080" ) )( 20 ), colorRampPalette( c( "#FF8080", "#FF0000" ) )( 30 ) ), limits=c( 0, 16 ), breaks=c( 0, 8, 16 ) ) +
	scale_size( range=c( 1, 7 ), breaks=c( 1, 4, 10 ), limits=c( 1, 10 ) ) +
	scale_x_continuous( breaks=seq( 1, length( GeneClass ) ), labels=seq( 1, length( GeneClass ) ), limits=c( 0, length( GeneClass ) + 1 ), expand=c( 0, 0 ) ) +
	scale_y_continuous( breaks=unlist( lapply( sort( unique( gg ) ), function(x) { return( median( match( names( which( gg == x ) ), go ) ) ); } ) ), 
		labels=sort( unique( gg ) ), limits=c( 0.5, length( gg ) + 0.5 ), expand=c( 0, 0 ) )+
	theme( legend.position="right", legend.justification="left", legend.background=element_rect( fill="#FFFFFF" ), legend.key=element_rect( fill="white" ), 
		aspect.ratio=78 / 34,
		plot.title=element_text( size=12 ), 
		panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_rect( colour="black", fill=NA, size=1.3 ), panel.background=element_blank(),
		axis.line.x=element_blank(), axis.line.y=element_blank(), 
		axis.text.x=element_text( colour="black", size=4.8, angle=90, hjust=1, vjust=0.5 ), axis.text.y=element_text( colour="black", hjust=1, vjust=0.5 ), 
		axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
		plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
	xlab( "" ) + ylab( "" ) + labs( title="" )
ggsave( file="tmp1.svg", plot=p4, width=10, height=5 )
# 	--> Figure 7g enrichment bubble 

output <- fit; 
output$ClusterName <- c( "A", "B", "C", "D", "E" )[ fit$Cluster ]; 
output$Enriched <- unlist( lapply( as.numeric( as.character( output$OR ) ), function(x) { return( ifelse( x > 1, "Y", "N" ) ); } ) ); 
output$SignificantEnriched <- unlist( apply( output, 1, function(x) { 
	
	if( as.numeric( as.character( x[ "OR" ] ) ) > 1 & as.numeric( as.character( x[ "BH" ] ) ) < 0.1 ) { return( "Y" ); }
	return( "N" ); 
	
	} ) ); 
write.table( output, file="tmp.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 
# 	--> Figure 7g enrichment bubble 

#Figure 7h 7i
library( "miceadds" );
library( "DescTools" )
library( "limma" )
library( "ggplot2" ); 
yTickFunc2 <- function(y) { 
	
	s <- 10^floor( log10( y ) ); 
	Exp_tick <- 0; 
	Exp_tick0 <- seq( 0, ceiling( y / s ) ) * s; 
	Exp_tick1 <- seq( 0, ceiling( y / s / 4 ) ) * s * 4; 
	Exp_tick2 <- seq( 0, ceiling( y / s / 5 ) ) * s * 5; 
	
	s <- 1; 
	Exp_tick6 <- seq( 0, ceiling( y / s ) ) * s; 
	Exp_tick7 <- seq( 0, ceiling( y / s / 4 ) ) * s * 4; 
	Exp_tick8 <- seq( 0, ceiling( y / s / 5 ) ) * s * 5; 
	
	Exp_tick <- list( Exp_tick0, Exp_tick1, Exp_tick2, Exp_tick6, Exp_tick7, Exp_tick8 )[ order( c( max( Exp_tick0 ) - y, max( Exp_tick1 ) - y, max( Exp_tick2 ) - y, max( Exp_tick6 ) - y, max( Exp_tick7 ) - y, max( Exp_tick8 ) - y ) ) ]; 
	return( unlist( Exp_tick[ which( unlist( lapply( Exp_tick, length ) ) %in% c( 3, 4, 5 ) )[ 1 ] ] ) ); 
	
}


sample_info <- read.csv( "sample_info.csv", header=TRUE, check.names=FALSE ); 

geneinfo <- read.table( "geneinfo.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 


# load ISG gene list 
ISG.gene <- as.vector( as.matrix( read.table( "ISG.gene.lst", header=FALSE, check.names=FALSE, sep="\t", quote="", comment.char="" ) ) );
print( length( ISG.gene ) ); 


t <- list(); 

t$ISG <- ISG.gene; 
geneSets <- t; 


# load expression from DepMap
load.Rdata( filename="DepMap_ExpFull.RData", "Exp" );
dfS <- cbind( Exp$dfS, data.frame( Tissue=unlist( lapply( as.character( Exp$dfS$CCLE_Name ), function(x) { 
	
	xx <- unlist( strsplit( x, "_" ) ); 
	return( paste( xx[ 2 : length( xx ) ], collapse="_" ) ); 
	
	} ) ) ) )

# manually correct sample list included for downstream analysis 
sample <- read.table( "sample.Exp.LZ.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
dfS$TissueCorrected <- as.character( sample[ match( as.character( dfS$DepMap_ID ), as.character( sample$DepMap_ID ) ), "L Zhang Cancer Types" ] )

dfS$Tissue <- as.character( dfS$TissueCorrected ); 
tab <- table( sort( as.character( dfS$Tissue ) ) ); 

# select cell lines to include 
ind <- which( as.character( dfS$TissueCorrected ) %in% setdiff( names( tab[ tab >= 3 ] ), c( "FIBROBLAST", "" ) ) ); 
df <- Exp$df[ ind, ]; 
dfS <- dfS[ ind, ]; 
tab <- table( sort( as.character( dfS$Tissue ) ) ); 

geneRow <- Exp$dfG; 


# calculate signatures of gene sets
t <- lapply( geneSets, function(gs) { 
	
	indR <- match( gs, as.character( geneRow$GeneName ) ); 
	indR <- indR[ ! is.na( indR ) ]; 
	Exp <- df[ , indR ]; 
	geneRow <- geneRow[ indR, ]; 
	colnames( Exp ) <- as.character( geneRow$GeneName ); 

	geneExpZ <- apply( Exp, 2, function(x) { 
		
		mad <- MeanAD( x ); 
		med <- median( x ); 
		return( ( x - med ) / ( 1.253314 * mad ) ); 
		
		} )
	t <- data.frame( t=unlist( apply( geneExpZ, 1, mean ) ) ); 
	return( t ); 
	
	} )
t <- do.call( cbind, t ); 
row.names( t ) <- row.names( df ); 
colnames( t ) <- names( geneSets ); 

sigs <- t; 


# association test 
t <- lapply( colnames( sigs ), function( gs ) { 
	
	sml <- sigs[ , gs ]; 

	ind <- which( ! is.na( sml ) ); 
	sml <- sml[ ind ]; 

	design <- model.matrix( ~ sml ); 
	
	fit <- lmFit( data.frame( t( df[ ind, ] ) ), design ); 
	fit <- eBayes( fit, robust=TRUE ); 
	
	tT <- subset( topTable( fit, adjust="BH", sort.by="p", number=ncol( df ) ), select=c( "adj.P.Val","P.Value","t","B","logFC" ) ); 

	t <- do.call( rbind, apply( df[ ind, row.names( tT ) ], 2, function(x) { 
		qt <- quantile( x, na.rm=TRUE ); 
		return( data.frame( qt3=qt[ "75%" ], qt1=qt[ "25%" ] ) ); 
	} ) )
	tT <- cbind( tT, t ); 
	
	tT$effect <- tT$logFC / sqrt( fit$s2.post[ row.names( tT ) ] ) * ( tT$qt3 - tT$qt1 ); 

	tT <- cbind( data.frame( GeneSet=gs, Gene=row.names( tT ) ), geneRow[ match( row.names( tT ), colnames( df ) ), ], tT ); 
	
	return( tT ); 
	
	} )
tT <- do.call( rbind, t ); 

t <- apply( tT, 1, function(x) { 
	
	if( as.character( x[ "GeneName" ] ) %in% geneSets[[ as.character( x[ "GeneSet" ] ) ]] ) { return( 1 ); 
	} else { return( 0 ); }
	
	} )
tT$Include <- unlist( t ); 


# annotate IAMs 

TE6 <- read.table( "CCLE.TE.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 

t <- subset( TE6, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` == 1 ); 
print( table( t$Group ) ); 
t$Detected <- unlist( lapply( t$Group, function(x) { return( ifelse( as.character( x ) %in% c( "Not-Expressed" ), "N", "Y" ) ); } ) ); 
t$Ubiquitous <- unlist( lapply( t$Group, function(x) { return( ifelse( as.character( x ) %in% c( "Expressed-In-All" ), "Y", "N" ) ); } ) ); 
t$Selective <- unlist( lapply( t$Group, function(x) { return( ifelse( as.character( x ) %in% c( "Bimodal-Like", "Lineage-Enriched", "Other", "Right-Skewed" ), "Y", "N" ) ); } ) ); 
t$Lineage <- unlist( lapply( t$Group, function(x) { return( ifelse( as.character( x ) %in% c( "Lineage-Enriched" ), "Y", "N" ) ); } ) ); 
write.table( t, file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 

MasterTableV31 <- read.table( "MasterTableV31.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
t <- subset( TE6, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` == 1 ); 
write.table( MasterTableV31[ match( t$GeneAcc, MasterTableV31$`Ensembl ID` ), ], file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 

general <- read.table( "general.txt", header=FALSE, check.names=FALSE, sep="\t", quote="", comment.char="" ); 
t <- subset( TE6, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` == 1 ); 
write.table( general[ match( t$GeneAcc, general$`V4` ), ], file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 


tmp <- as.character( subset( TE6, ! is.na( Surfaceome ) & `Immunological Accessory Molecules` == 1 )$GeneAcc )
tT$IAM <- unlist( lapply( as.character( tT$EnsemblID ), function(x) { return( ifelse( x %in% tmp, 1, 0 ) ); } ) )

tT$Group <- unlist( apply( tT, 1, function(x) { 
	return( as.character( subset( TE6, GeneAcc == as.character( x[ "EnsemblID" ] ) )$Group ) ); } ) )

tT$genetype <- as.character( geneinfo[ match( as.character( tT$EnsemblID ), as.character( geneinfo$ensg ) ), "genetype" ] )



FDRcut <- 10^-8; 


xmax <- max( abs( subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" )$logFC ), na.rm=TRUE ); 
xtick <- yTickFunc2( xmax ); 
ymax <- max( -log10( subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 )$P.Value ), na.rm=TRUE ); 
print( ymax ); 
print( head( subset( tT[ with( tT, order( abs( logFC ) ) ), ], ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 & IAM == 1 & adj.P.Val < FDRcut & logFC > 0 ) ) ); 
tT_red <- subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 & IAM == 1 & adj.P.Val < FDRcut & logFC > 0.25 ); 
p_Exp <- ggplot(  ) + 
	geom_point( data=subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 & IAM == 0 ), 
		aes( logFC, -log10( P.Value ) ), shape=16, col="#9E9E9E", alpha=0.8, size=1.8 ) +
	geom_point( data=subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 & IAM == 1 & ! Gene %in% tT_red$Gene ), 
		aes( logFC, -log10( P.Value ) ), shape=16, col="#0000FF", alpha=0.8, size=1.8 ) +
	geom_point( data=tT_red, 
		aes( logFC, -log10( P.Value ) ), shape=16, col="#FF0000", alpha=0.8, size=1.8 ) +
	scale_x_continuous( breaks=c( -4, -2, 0, 2, 4 ), limits=c( -max( xtick ), max( xtick ) ) ) +
	scale_y_continuous( breaks=seq( 0, 250, 50 ), limits=c( 0, 250 ) ) +
	scale_size_continuous( name="", range=c( 0.5, 3 ) ) +
	geom_hline( yintercept=-log10( min( subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & ! is.na( adj.P.Val ) & adj.P.Val >= FDRcut )$P.Value ) ), lty=2, col="#BFBFBF" ) +
	geom_vline( xintercept=0.25, lty=2, col="#BFBFBF" ) +
	geom_vline( xintercept=-0.25, lty=2, col="#BFBFBF" ) +
	xlab( "" ) + ylab( "" ) + ggtitle( "Protein coding genes (Expressed)" ) +
	theme( legend.position="none", aspect.ratio=1, legend.key=element_blank(), 
		axis.line.x=element_line( size=0.7, colour="black" ), axis.line.y=element_line( size=0.7, colour="black" ), 
		panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), 
		axis.ticks.x=element_line( colour="black", size=0.6, linetype=1 ), axis.ticks.y=element_line( colour="black", size=0.6, linetype=1 ), axis.ticks.length=unit( 0.14, "cm" ) ); 
# pdf( "tmp2.pdf", width=10, height=5 ); 
# print( p_Exp ); 
# dev.off()
ggsave( file="tmp2.svg", plot=p_Exp, width=10, height=5 )
# 	--> Figure 7H volcano

t1 <- subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 & IAM == 0 ); 
t1$GeneDefine <- "Other"
t2 <- subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 & IAM == 1 & ! Gene %in% tT_red$Gene ); 
t2$GeneDefine <- "mIAM, not IFN stimulating"
t3 <- tT_red; 
t3$GeneDefine <- "mIAM, IFN stimulating"; 
output <- rbind( t1, t2, t3 ); 
output$Y <- -log10( output$P.Value ); 
output$mIAM <- c( "N", "Y" )[ match( output$IAM, c( 0, 1 ) ) ]; 
write.table( output, file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 
# 	--> Figure 7H volcano


# calculate enrichment, IAM / non-IAM ~ red (exclude genes included in ISG signature)

dfT <- subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 ); 
dfT$Gene <- factor( c( "Other PcG", "IAM" )[ match( as.character( dfT$IAM ), c( "0", "1" ) ) ], levels=c( "Other PcG", "IAM" ) ); 

FDRcut <- 10^-8; 
GI <- as.character( subset( dfT, adj.P.Val < FDRcut & logFC > 0 )$EnsemblID )
dfT$GI <- factor( as.character( dfT$EnsemblID ) %in% GI, levels=c( FALSE, TRUE ) )

tab <- table( dfT$Gene, dfT$GI ); 
print( tab ); 
write.table( as.data.frame( tab ), file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 
fit <- fisher.test( tab ); 
t1 <- data.frame( GeneClass="IAM", Group="Group I", OR=as.numeric( as.character( fit$estimate ) ), pvalue=as.numeric( as.character( fit$p.value ) ) ); 
print( t1 ); 
# IAM ~ IFN stimulating, OR=1.535, pvalue=0.000186

Pct <- do.call( rbind, lapply( c( "Other PcG", "IAM" ), function(GRPi) { 
	
	t <- subset( dfT, Gene == GRPi ); 
	t$GI <- factor( as.character( t$GI ), levels=c( FALSE, TRUE ) ); 
	r <- cbind( data.frame( Gene=GRPi, N=nrow( t ) ), as.data.frame( table( t$GI ) / nrow( t ) ) ); 
	colnames( r )[ 3:4 ] <- c( "GI", "Pct" ); 
	return( r ); 
	
	} ) )
Pct <- do.call( rbind, lapply( unique( as.character( Pct$Gene ) ), function(x) { 
	
	cat( x, "; " ); 
	t <- subset( Pct, Gene == x ); 
	t <- do.call( rbind, lapply( seq( 1, nrow( t ) ), function(y) { 
		
		ymin <- ifelse( y == 1, 0, sum( t[ 1:(y-1), ]$Pct ) ); 
		ymax <- ymin + t[ y, ]$Pct; 
		return( cbind( t[ y, ], data.frame( ymin=ymin, ymax=ymax ) ) ); 
		
		} ) ); 
	return( t ); 
	
	} ) )
p3_mosaic <- ggplot( ) +
	geom_rect( data=subset( Pct, Gene=="IAM" ), aes( fill=as.character( GI ), xmin=1 - 0.4, xmax=1 + 0.4, ymin=ymin, ymax=ymax ) ) +
	geom_rect( data=subset( Pct, Gene=="Other PcG" ), aes( fill=as.character( GI ), xmin=0 - 0.4, xmax=0 + 0.4, ymin=ymin, ymax=ymax ) ) +
	scale_x_continuous( limits=c( -1, 2 ), breaks=NULL, labels=NULL ) +
	coord_polar( theta="y" ) +
	xlab( "" ) + ylab( "" ) + ggtitle( "" ) +
	scale_fill_manual( values=c( "TRUE"="#FFA500", "FALSE"="#D9D9D9" ) ) + 
	theme( legend.position="none", aspect.ratio=1, 
		panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
		panel.border=element_blank(), panel.background=element_blank(), 
		axis.text.x=element_blank(), axis.text.y=element_blank(), 
		axis.line.x=element_line( size=0.8, colour="black" ), axis.line.y=element_line( size=0.8, colour="black" ), 
		axis.ticks.x=element_line( colour="black", size=0.5, linetype=1 ), axis.ticks.y=element_line( colour="black", size=0.5, linetype=1 ), axis.ticks.length=unit( 0.12, "cm" ), 
		plot.margin=unit( c( 0, 0, 0, 1 ), "cm" ) )
ggsave( file="tmp.svg", plot=p3_mosaic, width=3, height=5 )
# 	--> Figure 7H donut


# calculate enrichment, within IAM, category ~ red 

dfT <- subset( tT, ! is.na( genetype ) & genetype == "protein_coding" & Group != "Not-Expressed" & Include == 0 & IAM == 1 ); 

GI <- as.character( subset( dfT, adj.P.Val < FDRcut & logFC > 0 )$EnsemblID )
dfT$GI <- factor( as.character( dfT$EnsemblID ) %in% GI, levels=c( FALSE, TRUE ) )

dfAnnoImmune <- read.table( "immune_surface.ImmuneV29.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 

# GeneClass <- c( "Adhesion Molecules", "Fc Receptors", "Tetraspanins", "Cytokine", "Cytokine Receptors", "Chemokine", "Chemokine Receptors", "Co-stimulator/Co-inhibitor", "Membrane-bound pattern recognition receptors" ); 
GeneClass <- c( "Adhesion Molecules", "Fc Receptors", "Tetraspanins", "Cytokine/Cytokine Receptors", "Chemokine/Chemokine Receptors", "Co-stimulator/Co-inhibitor", "Membrane-bound pattern recognition receptors" ); 
t <- lapply( seq( 1, length( GeneClass ) ), function(i) { 
	
	geneClass <- GeneClass[ i ]; 
	
	if( geneClass == "Cytokine/Cytokine Receptors" ) { 
		genes <- unique( dfAnnoImmune[ as.character( dfAnnoImmune[ , "Cytokine" ] ) != "", "Ensembl ID" ], dfAnnoImmune[ as.character( dfAnnoImmune[ , "Cytokine Receptors" ] ) != "", "Ensembl ID" ] ); 
	} else if( geneClass == "Chemokine/Chemokine Receptors" ) { 
		genes <- unique( dfAnnoImmune[ as.character( dfAnnoImmune[ , "Chemokine" ] ) != "", "Ensembl ID" ], dfAnnoImmune[ as.character( dfAnnoImmune[ , "Chemokine Receptors" ] ) != "", "Ensembl ID" ] ); 
	} else { 
		genes <- dfAnnoImmune[ as.character( dfAnnoImmune[ , geneClass ] ) != "", "Ensembl ID" ]; 
	}
	
	tab <- table( dfT$GI, dfT$EnsemblID %in% genes ); 
	fit <- fisher.test( tab ); 
	t1 <- data.frame( GeneClass=geneClass, Group="Group I", OR=as.numeric( as.character( fit$estimate ) ), pvalue=as.numeric( as.character( fit$p.value ) ) ); 
	print( t1 ); 
	
	return( t1 ); 
	
	} )
fit <- do.call( rbind, t ); 

fit$BH <- p.adjust( fit$pvalue, method="BH" ); 
fit$variable <- as.numeric( factor( as.character( fit$GeneClass ), levels=rev( GeneClass ) ) );
fit$significance <- unlist( lapply( fit$pvalue, function(x) { return( ifelse( x <= 1e-06, 1e-06, x ) ); } ) )
print( sort( round( log2( fit$OR ), 2 ) ) ); 
print( fit ); 
write.table( fit, file="tmpS1.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 
# 	--> Figure 7I bar

p3_1 <- ggplot(  ) +
	geom_bar( data=fit, aes( variable, -log10( significance ) * ifelse( OR > 1, 1, -1 ), 
		fill=ifelse( abs( log2( OR ) ) > 1.6, 1.6 * sign( log2( OR ) ), log2( OR ) ) ), color=NA, stat="identity", width=0.8 ) +
	geom_segment( aes( x=0.2, xend=( nrow( fit ) + 0.8 ), y=0, yend=0 ), size=0.5, col="#000000" ) +
	scale_fill_gradient2( guide=guide_colourbar( barwidth=0.5, barheight=3 ), name="log2(OR)", low="#B07405", mid="white", high="#662d91", midpoint=0, space="rgb", na.value="grey80", 
		limits=c( -1.6, 1.6 ), breaks=c( -1.6, -0.8, 0, 0.8, 1.6 ) ) +
	xlab( "" ) + ylab( "" ) + labs( title="" ) +
	scale_x_continuous( breaks=c( seq( 1, length( GeneClass ) ) ), labels=unlist( lapply( rev( GeneClass ), function(x) { return( substr( x, 1, 28 ) ); } ) ), expand=c( 0, 0 ) ) +
	scale_y_continuous( expand=c( 0, 0 ) ) +
	coord_flip( xlim=c( 0.2, ( nrow( fit ) + 0.8 ) ), ylim=c( -6, 6 ) ) +
	theme( legend.position="right", aspect.ratio=5 / 3, 
		axis.line.x=element_line( size=0.5, colour="black" ), axis.line.y=element_blank(), 
		panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
		panel.border=element_blank(), 
		axis.ticks.x=element_line( colour="black", size=0.4, linetype=1 ), axis.ticks.y=element_blank(), axis.ticks.length=unit( 0.09, "cm" ) )
ggsave( file="tmp.svg", plot=p3_1, width=4, height=5 )
# 	--> Figure 7I bar



