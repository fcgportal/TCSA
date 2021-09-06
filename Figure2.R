#Figure 2d
library( "ggplot2" ); 

expMedian <- read.table( "Fig2_expMedian.txt", header=TRUE, check.names=FALSE, row.names=1, sep="\t", quote="", comment.char="" )

tlst <- colnames( expMedian ); 
tlst0 <- tlst[ grepl( "^NormalGTEX", tlst ) ]
tlst1 <- tlst[ grepl( "^TCGA-", tlst ) ]; 


# load Column DO [Subcellular in Surfaceome paper] of master table V39
SubLoc <- read.table( "Fig2_Surfaceome.txt", header=TRUE, check.names=FALSE, sep="\t", quote="", comment.char="#" ); 
print( dim( SubLoc ) ); # 60,498
# print( showd( SubLoc ) ); 
print( table( SubLoc$`Subcellular in Surfaceome paper` ) ); 


# all subcellular locations
gorder <- c( "Surface", "Endoplasmic_reticulum", "Nucleus", "Golgi_apparatus", "Cytoplasm", "Peroxisome", "Lysosome/Vacuole", "Mitochondrion" ); 
print( table( gorder %in% names( table( SubLoc$`Subcellular in Surfaceome paper` ) ) )  ); 


# for each gene, calculate number of tissues (TCGA) in which it is expressed (median FPKM > 1)
SubLoc$SubcellularLocalization <- as.character( SubLoc$`Subcellular in Surfaceome paper` ); 
SubLoc$ENSG <- as.character( SubLoc$`Ensembl ID` ); 
SubLoc <- subset( SubLoc, SubcellularLocalization != "" ); 
print( table( as.character( SubLoc$ENSG ) %in% row.names( expMedian ) ) )

SubLoc$NumTissue <- apply( expMedian[ as.character( SubLoc$ENSG ), tlst1 ], 1, function(x) { return( length( which( x > 1 ) ) ); } )
SubLoc$NumTissueBin <- unlist( lapply( SubLoc$NumTissue, function(x) { return( ifelse( x >= 26, "U", ifelse( x <= 6, "S", "O" ) ) ); } ) ); 

dfp <- do.call( rbind, lapply( gorder, function(g) { 
	
	t1 <- as.character( subset( SubLoc, SubcellularLocalization == g )$NumTissueBin ); 
	t <- cbind( data.frame( group=g ), as.data.frame( table( factor( t1, levels=c( "S", "O", "U" ) ) ) ), as.data.frame( table( factor( t1, levels=c( "S", "O", "U" ) ) ) / length( t1 ) ) ); 
	colnames( t ) <- c( "group", "Var1", "GeneNum", "ExpGroup", "Freq" ); 
	return( t ); 
	
	} ) )

# reorder subcellular locations (order by percentage of specific genes)
print( gorder ); 
gorder <- as.character( subset( dfp[ with( dfp, rev( order( Freq ) ) ), ], Var1 == "S" )$group ); 
print( gorder ); 
cat( paste( gorder, collapse="\", \"" ) ); 
# c( "Surface", "Lysosome/Vacuole", "Endoplasmic_reticulum", "Nucleus", "Golgi_apparatus", "Mitochondrion", "Cytoplasm", "Peroxisome" )

print( max( dfp$Freq ) ); 

p2 <- ggplot( data=subset( dfp, group %in% gorder ), aes( x=Var1, y=Freq, fill=factor( group, levels=gorder ) ) ) +
	geom_bar( width=0.7, stat="identity", color="#7F7F7F", position=position_dodge() )+
	scale_fill_manual( name="", values=c( "Surface"="#FF0000", "Lysosome/Vacuole"="#0040FF", "Endoplasmic_reticulum"="#0095FF", "Nucleus"="#00EAFF", 
		"Golgi_apparatus"="#00FF6A", "Mitochondrion"="#AAFF00", "Cytoplasm"="#FFFF00", "Peroxisome"="#FFD500" ) ) +
	scale_x_discrete( breaks=c( "S", "O", "U" ), labels=c( "0-6", "7-25", "26-33" ) ) + 
	scale_y_continuous( breaks=c( 0, 0.2, 0.4, 0.6, 0.8 ), limits=c( 0, 0.8 ) ) +
	labs( title="" ) + xlab( "# of tissues" ) + ylab( "% of genes (FPKM > 1)" ) + 
	theme( legend.position="right", aspect.ratio=1 / 1.9, 
		axis.line.x=element_line( colour="black", size=1.2, linetype=1 ), axis.line.y=element_line( colour="black", size=1.2, linetype=1 ), 
		panel.border=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
		axis.ticks.x=element_line( colour="black", size=1, linetype=1 ), axis.ticks.y=element_line( colour="black", size=1, linetype=1 ), axis.ticks.length=unit( 0.2, "cm" ), 
		plot.margin=unit( c( 0, 0, 0, 0 ), "cm" ) )
svg( filename="Figure2D.svg", width=8.4, height=4.2 ); print( p2 ); dev.off(); 

write.table( subset( dfp, group %in% gorder ), file="Figure2D.txt", append=FALSE, quote=F, row.names=FALSE, col.names=TRUE, sep="\t" ); 


#Figure 2g
library( "tsne" )
library( "ggplot2" )

tsne <- read.table("Fig2_AllSamples.tSNE.txt", header=TRUE, check.names=FALSE, row.names=1, sep="\t", quote="", comment.char="" );
colnames( tsne ) <- c( "tSNE1", "tSNE2" )
row.names( tsne ) <- as.character( unlist( RPKM[ 2, col.index ] ) )
tsne[ , "cancertype" ] <- as.character( unlist( RPKM[ 1, col.index ] ) )
p <- ggplot(  ) + 
		# geom_text( data=tsne, aes( x=tSNE1, y=tSNE2, label=as.character( cancertype ), color=as.character( cancertype ) ), size=1.5, alpha=0.5 ) + 
		geom_point( data=tsne, aes( x=tSNE1, y=tSNE2, label=as.character( cancertype ), color=as.character( cancertype ) ), size=1.5, alpha=0.5 ) +
		scale_color_manual( values=c("Tumor-Adrenal_Gland-ACC"="#B03A2E","Tumor-Bladder-BLCA"="#873600","Tumor-Breast-BRCA"="#DEF2C7","Tumor-Cervix-CESC"="#74B02E","Tumor-Bile_Duct-CHOL"="#C0392B","Tumor-Colon-COAD"="#FDEDEC","Tumor-Blood-DLBC"="#76CCD5","Tumor-Esophagus-ESCA"="#633974","Tumor-Brain-GBM"="#35AFBB","Tumor-Head_Neck-HNSC"="#76448A","Tumor-Kidney-KICH"="#E59866","Tumor-Kidney-KIRC"="#FAE5D3","Tumor-Kidney-KIRP"="#F5CBA7","Tumor-Blood-LAML"="#CBF1F4","Tumor-Brain-LGG"="#2CD8E9","Tumor-Liver-LIHC"="#E74C3C","Tumor-Lung-LUAD"="#D2B4DE","Tumor-Lung-LUSC"="#9B59B6","Tumor-Body_Cavities-MESO"="#E8DAEF","Tumor-Ovary-OV"="#CBE4AE","Tumor-Pancreas-PAAD"="#E6B0AA","Tumor-Paraganglia-PCPG"="#641E16","Tumor-Prostate-PRAD"="#D35400","Tumor-Colon-READ"="#F2D7D5","Tumor-Soft_Tissue-SARC"="#415A24","Tumor-Skin-SKCM"="#075B63","Tumor-Stomach-STAD"="#4A235A","Tumor-Testis-TGCT"="#D4E6F1","Tumor-Thyroid-THCA"="#E8F8F5","Tumor-Thymus-THYM"="#F5EEF8","Tumor-Uterus-UCEC"="#A1C27C","Tumor-Uterus-UCS"="#808000","Tumor-Eye-UVM"="#048693","Adjacent-Bladder-BLCA"="#DCDCDC","Adjacent-Breast-BRCA"="#DCDCDC","Adjacent-Cervix-CESC"="#DCDCDC","Adjacent-Bile_Duct-CHOL"="#DCDCDC","Adjacent-Colon-COAD"="#DCDCDC","Adjacent-Esophagus-ESCA"="#DCDCDC","Adjacent-Brain-GBM"="#DCDCDC","Adjacent-Head_Neck-HNSC"="#DCDCDC","Adjacent-Kidney-KICH"="#DCDCDC","Adjacent-Kidney-KIRC"="#DCDCDC","Adjacent-Kidney-KIRP"="#DCDCDC","Adjacent-Liver-LIHC"="#DCDCDC","Adjacent-Lung-LUAD"="#DCDCDC","Adjacent-Lung-LUSC"="#DCDCDC","Adjacent-Pancreas-PAAD"="#DCDCDC","Adjacent-Paraganglia-PCPG"="#DCDCDC","Adjacent-Prostate-PRAD"="#DCDCDC","Adjacent-Colon-READ"="#DCDCDC","Adjacent-Soft_Tissue-SARC"="#DCDCDC","Adjacent-Skin-SKCM"="#DCDCDC","Adjacent-Stomach-STAD"="#DCDCDC","Adjacent-Thyroid-THCA"="#DCDCDC","Adjacent-Thymus-THYM"="#DCDCDC","Adjacent-Uterus-UCEC"="#DCDCDC","Normal-Adipose_Tissue"="#DCDCDC","Normal-Muscle"="#DCDCDC","Normal-Blood_Vessel"="#DCDCDC","Normal-Heart"="#DCDCDC","Normal-Ovary"="#DCDCDC","Normal-Uterus"="#DCDCDC","Normal-Breast"="#DCDCDC","Normal-Salivary_Gland"="#DCDCDC","Normal-Brain"="#DCDCDC","Normal-Adrenal_Gland"="#DCDCDC","Normal-Thyroid"="#DCDCDC","Normal-Lung"="#DCDCDC","Normal-Pancreas"="#DCDCDC","Normal-Esophagus"="#DCDCDC","Normal-Stomach"="#DCDCDC","Normal-Skin"="#DCDCDC","Normal-Colon"="#DCDCDC","Normal-Small_Intestine"="#DCDCDC","Normal-Prostate"="#DCDCDC","Normal-Testis"="#DCDCDC","Normal-Nerve"="#DCDCDC","Normal-Spleen"="#DCDCDC","Normal-Pituitary"="#DCDCDC","Normal-Blood"="#DCDCDC","Normal-Vagina"="#DCDCDC","Normal-Liver"="#DCDCDC","Normal-Kidney"="#DCDCDC","Normal-Bladder"="#DCDCDC","Normal-Fallopian_Tube"="#DCDCDC","Normal-Cervix"="#DCDCDC","" ) ) +
		xlab( "tSNE1" ) + 
		ylab( "tSNE2" ) +
		ggtitle( "" ) +
		theme( legend.position="bottom", aspect.ratio=1, 
			panel.background=element_blank(), 
			panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
			panel.border=element_rect( colour="black", fill=NA, size=3.5 ), 
			axis.line.x=element_blank(), axis.line.y=element_blank(), 
			axis.text.x=element_blank(), axis.text.y=element_blank(), 
			axis.ticks.x=element_line( colour="black", size=1, linetype=1 ), axis.ticks.length=unit( 0.25, "cm" ), 
			axis.ticks.y=element_line( colour="black", size=1, linetype=1 ), 
			plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
		guides( size=FALSE )
svg( "Fig2_AllSamples.tSNE.svg", width=10, height=15 ); print( p ); dev.off()

wss <- sapply( 17:51, function(k){ 
	kmeans( tsne[ , c( "tSNE1", "tSNE2" ) ], k, nstart=50, iter.max=15 )$tot.withinss } )
svg( "Fig2_AllSamples.tSNE(kmeans).svg" ); plot( seq( 17, 51 ), wss ); dev.off()

k <- 33
tsne.kmeans <- kmeans( tsne[ , c( "tSNE1", "tSNE2" ) ], k, nstart=50, iter.max=15 )
fit <- chisq.test( table( tsne$cancertype, tsne.kmeans$cluster ) )
print( pchisq( fit$statistic, df=fit$parameter, lower.tail=FALSE, log.p=TRUE ) )

tsne <- read.table("Fig2_Prostate.tSNE.txt", header=TRUE, check.names=FALSE, row.names=1, sep="\t", quote="", comment.char="" );
colnames( tsne ) <- c( "tSNE1", "tSNE2" )
row.names( tsne ) <- as.character( unlist( RPKM[ 2, col.index ] ) )
tsne[ , "cancertype" ] <- as.character( unlist( RPKM[ 1, col.index ] ) )
p <- ggplot(  ) + 
		# geom_text( data=tsne, aes( x=tSNE1, y=tSNE2, label=as.character( cancertype ), color=as.character( cancertype ) ), size=1.5, alpha=0.5 ) + 
		geom_point( data=tsne, aes( x=tSNE1, y=tSNE2, label=as.character( cancertype ), color=as.character( cancertype ) ), size=5, alpha=0.5 ) +
		scale_color_manual( values=c("Tumor-Adrenal_Gland-ACC"="#DCDCDC","Tumor-Bladder-BLCA"="#DCDCDC","Tumor-Breast-BRCA"="#DCDCDC","Tumor-Cervix-CESC"="#DCDCDC","Tumor-Bile_Duct-CHOL"="#DCDCDC","Tumor-Colon-COAD"="#DCDCDC","Tumor-Blood-DLBC"="#DCDCDC","Tumor-Esophagus-ESCA"="#DCDCDC","Tumor-Brain-GBM"="#DCDCDC","Tumor-Head_Neck-HNSC"="#DCDCDC","Tumor-Kidney-KICH"="#DCDCDC","Tumor-Kidney-KIRC"="#DCDCDC","Tumor-Kidney-KIRP"="#DCDCDC","Tumor-Blood-LAML"="#DCDCDC","Tumor-Brain-LGG"="#DCDCDC","Tumor-Liver-LIHC"="#DCDCDC","Tumor-Lung-LUAD"="#DCDCDC","Tumor-Lung-LUSC"="#DCDCDC","Tumor-Body_Cavities-MESO"="#DCDCDC","Tumor-Ovary-OV"="#DCDCDC","Tumor-Pancreas-PAAD"="#DCDCDC","Tumor-Paraganglia-PCPG"="#DCDCDC","Tumor-Prostate-PRAD"="#D35400","Tumor-Colon-READ"="#DCDCDC","Tumor-Soft_Tissue-SARC"="#DCDCDC","Tumor-Skin-SKCM"="#DCDCDC","Tumor-Stomach-STAD"="#DCDCDC","Tumor-Testis-TGCT"="#DCDCDC","Tumor-Thyroid-THCA"="#DCDCDC","Tumor-Thymus-THYM"="#DCDCDC","Tumor-Uterus-UCEC"="#DCDCDC","Tumor-Uterus-UCS"="#DCDCDC","Tumor-Eye-UVM"="#DCDCDC","Adjacent-Bladder-BLCA"="#DCDCDC","Adjacent-Breast-BRCA"="#DCDCDC","Adjacent-Cervix-CESC"="#DCDCDC","Adjacent-Bile_Duct-CHOL"="#DCDCDC","Adjacent-Colon-COAD"="#DCDCDC","Adjacent-Esophagus-ESCA"="#DCDCDC","Adjacent-Brain-GBM"="#DCDCDC","Adjacent-Head_Neck-HNSC"="#DCDCDC","Adjacent-Kidney-KICH"="#DCDCDC","Adjacent-Kidney-KIRC"="#DCDCDC","Adjacent-Kidney-KIRP"="#DCDCDC","Adjacent-Liver-LIHC"="#DCDCDC","Adjacent-Lung-LUAD"="#DCDCDC","Adjacent-Lung-LUSC"="#DCDCDC","Adjacent-Pancreas-PAAD"="#DCDCDC","Adjacent-Paraganglia-PCPG"="#DCDCDC","Adjacent-Prostate-PRAD"="#048693","Adjacent-Colon-READ"="#DCDCDC","Adjacent-Soft_Tissue-SARC"="#DCDCDC","Adjacent-Skin-SKCM"="#DCDCDC","Adjacent-Stomach-STAD"="#DCDCDC","Adjacent-Thyroid-THCA"="#DCDCDC","Adjacent-Thymus-THYM"="#DCDCDC","Adjacent-Uterus-UCEC"="#DCDCDC","Normal-Adipose_Tissue"="#DCDCDC","Normal-Muscle"="#DCDCDC","Normal-Blood_Vessel"="#DCDCDC","Normal-Heart"="#DCDCDC","Normal-Ovary"="#DCDCDC","Normal-Uterus"="#DCDCDC","Normal-Breast"="#DCDCDC","Normal-Salivary_Gland"="#DCDCDC","Normal-Brain"="#DCDCDC","Normal-Adrenal_Gland"="#DCDCDC","Normal-Thyroid"="#DCDCDC","Normal-Lung"="#DCDCDC","Normal-Pancreas"="#DCDCDC","Normal-Esophagus"="#DCDCDC","Normal-Stomach"="#DCDCDC","Normal-Skin"="#DCDCDC","Normal-Colon"="#DCDCDC","Normal-Small_Intestine"="#DCDCDC","Normal-Prostate"="#A1C27C","Normal-Testis"="#DCDCDC","Normal-Nerve"="#DCDCDC","Normal-Spleen"="#DCDCDC","Normal-Pituitary"="#DCDCDC","Normal-Blood"="#DCDCDC","Normal-Vagina"="#DCDCDC","Normal-Liver"="#DCDCDC","Normal-Kidney"="#DCDCDC","Normal-Bladder"="#DCDCDC","Normal-Fallopian_Tube"="#DCDCDC","Normal-Cervix"="#DCDCDC","" ) ) +
		xlab( "tSNE1" ) + 
		ylab( "tSNE2" ) +
		ggtitle( "" ) +
		theme( legend.position="bottom", aspect.ratio=1, 
			panel.background=element_blank(), 
			panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
			panel.border=element_rect( colour="black", fill=NA, size=3.5 ), 
			axis.line.x=element_blank(), axis.line.y=element_blank(), 
			axis.text.x=element_blank(), axis.text.y=element_blank(), 
			axis.ticks.x=element_line( colour="black", size=1, linetype=1 ), axis.ticks.length=unit( 0.25, "cm" ), 
			axis.ticks.y=element_line( colour="black", size=1, linetype=1 ), 
			plot.margin=unit(c( 0, 0, 0, 0 ), "cm" ) ) +
		guides( size=FALSE )
svg( "Fig2_Prostate.tSNE.svg", width=10, height=15 ); print( p ); dev.off()

wss <- sapply( 17:51, function(k){ 
	kmeans( tsne[ , c( "tSNE1", "tSNE2" ) ], k, nstart=50, iter.max=15 )$tot.withinss } )
svg( "Fig2_Prostate.tSNE(kmeans).svg" ); plot( seq( 17, 51 ), wss ); dev.off()

k <- 33
tsne.kmeans <- kmeans( tsne[ , c( "tSNE1", "tSNE2" ) ], k, nstart=50, iter.max=15 )
fit <- chisq.test( table( tsne$cancertype, tsne.kmeans$cluster ) )
print( pchisq( fit$statistic, df=fit$parameter, lower.tail=FALSE, log.p=TRUE ) )

#Figure 2k
library( "reshape2" )
library( "ggplot2" )
library( "grid" )

Percentage <- read.table( "Fig2_NES.txt", sep="\t", quote="", header=TRUE, check.names=FALSE )
Vote <- read.table( "Fig2_logp-val.txt", sep="\t", quote="", header=TRUE, check.names=FALSE )

gene.lst <- as.character( Percentage[ , "gene" ] )
gene.shown.lst <- make.names( gene.lst, unique=TRUE )
Percentage <- cbind( data.frame( gene=gene.shown.lst ), Percentage[ , 3 : ncol( Percentage ) ] )
Vote <- cbind( data.frame( gene=gene.shown.lst ), Vote[ , 3 : ncol( Vote ) ] )

study.lst <- colnames( Percentage )[ 2 : ncol( Percentage ) ]

df <- cbind( melt( Percentage, id="gene" ), melt( Vote, id="gene" )[ , "value", drop=FALSE ] )
colnames( df ) <- c( "gene", "study", "percentage", "vote" )


df[ , "x" ] <- match( df[ , "study" ], study.lst )
df[ , "y" ] <- match( df[ , "gene" ], rev( gene.shown.lst ) )
df[ , "color" ] <- ifelse( df$vote<1.301029996, "gainsboro", "red" )
df[ , "alpha" ] <- ifelse( df$vote<1.301029996, 1, df$vote / 3.000434077 )


svg( "Fig2_Enrichment.svg", width=12, height=10 )
layout( matrix( c( 1, 1, 1, 2, 3, 4 ), 3, 2, byrow=FALSE ) )

symbols( x=df$x, y=df$y, circles=( df$percentage ), inches=0.5 / 2, ann=F, bg=alpha( df$color, df$alpha ), fg="white", axes=FALSE )
axis( 1, at=c( 1:length( study.lst ) ), labels=study.lst, tick=FALSE, pos=0, hadj=1, las=3 ) # size = 1 means below
axis( 2, at=c( 1:length( gene.shown.lst ) ), labels=rev( gene.shown.lst ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 

size_legend <- signif( quantile( as.numeric( df$percentage ) ), 1 )
par( pin=c( 7 / 2, 1 / 2 ) )
symbols( x=factor( seq( 1, 6 ) ), y=rep( 1, 6 ), circles=( c( size_legend[ 1 ], 0.5, 0.75, 1, 1.5, size_legend[ 5 ] ) ), inches=0.5 / 2, ann=F, bg="red", fg="white", axes=FALSE )
axis( 2, at=1, labels="NES", tick=FALSE, pos=0, hadj=1, las=1 ) # size = 2 means left 
axis( 1, at=c( 1:6 ), labels=paste( as.character( c( size_legend[ 1 ], 0.5, 0.75, 1, 1.5, size_legend[ 5 ] ) ), "", sep="" ), tick=FALSE, pos=0, hadj=1, las=1 ) # size = 1 means below


vote_legend <- c( 1.301029996, 2, 3, 4)
par( pin=c( 4 / 2, 1 / 2 ) )
symbols( x=factor( seq( 1, 4 ) ), y=rep( 1, 4 ), circles=rep( 15, 4 ), inches=0.25 / 2, ann=F, bg=alpha( rep( "red", 4 ), ( vote_legend ) / 4 ), fg="white", axes=FALSE )
axis( 2, at=1, labels="p-val", tick=FALSE, pos=0, hadj=1, las=1 ) # vote = 2 means left 
axis( 1, at=c( 1:3 ), paste( as.character( c( 0.05, 0.01, 0.001 ) ), "", sep="" ), tick=FALSE, pos=0, hadj=1, las=1 ) 

dev.off()


