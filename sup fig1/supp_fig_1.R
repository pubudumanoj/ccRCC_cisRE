## wget("https://egg2.wustl.edu/roadmap/data/byFileType/peaks/unconsolidated/broadPeak/BI.Adult_Kidney.H3K4me1.153.broadPeak.gz")

## ------------------------------------------------------------------------


#first need to convert them to grc38. they are mapped to hg19

	#module load gcc/5.4.0 r-bundle-bioconductor/3.4
	
	#R
	
	#first need to convert HG19 to HG38 over chain
	
library(dplyr)
library(data.table)
library(stringi)
library(GenomicRanges)
library(rtracklayer)
hg19togrc38 <- function (fileN=NULL){
#load chain file
chain <- import.chain("../hg19ToHg38.over.chain")

#load subject file
a <- fread(file=fileN)

#select columns
a <- a[,c(1:3)]
colnames(a) <- c("chr","start","end")

df<- makeGRangesFromDataFrame(a, keep.extra.columns=TRUE)

#convert to hg38    
tx_hg19 <- liftOver(df, chain) %>% as.data.frame
	
	 tx_hg19 <- select(tx_hg19, c("seqnames","start","end"))
	 return(tx_hg19)
}


	 write.table(hg19togrc38(tx_hg19, fileN="BI.Adult_Kidney.H3K4me1.153.broadPeak"), file="BI.Adult_Kidney.H3K4me1.153.grc38.broadPeak", sep="\t", col.names=F, row.names=F, quote=F)
	 
	 write.table(hg19togrc38(tx_hg19, fileN="BI.Adult_Kidney.H3K27ac.153.broadPeak"), file="BI.Adult_Kidney.H3K27ac.153.grc38.broadPeak", sep="\t", col.names=F, row.names=F, quote=F)
	 
	 rite.table(hg19togrc38(tx_hg19, fileN="BI.Adult_Kidney.H3K4me3.153.broadPeak"), file="BI.Adult_Kidney.H3K4me3.153.grc38.broadPeak", sep="\t", col.names=F, row.names=F, quote=F)
	

#	q()
#	n
	

## 
## ------------------------------------------------------------------------
library(GenomicRanges)
library(rtracklayer)


## ------------------------------------------------------------------------

chi_overlap_barplot <- function(df=NULL, mark=NULL, random_ref=F){
if(random_ref==F){
  ref_h3k27ac <- fread(paste0("C:/Users/pubud/Documents/R_analysis/chip/roadmap peaks/BI.Adult_Kidney.",mark,".153.grc38.broadPeak"))
} else{
  ref_h3k27ac <- fread(paste0("C:/Users/pubud/Documents/R_analysis/chip/roadmap peaks/BI.Adult_Kidney.",mark,".153.grc38.random_seed1600.broadPeak"))
}
  
ref_h3k27ac <- select(ref_h3k27ac, c("V1","V2","V3"))
colnames(ref_h3k27ac) <- c("chr", "start", "end" )
	chr_list <- c(paste0(rep("chr", 22),1:22), "chrX","chrY")


ref_h3k27ac <- ref_h3k27ac[ref_h3k27ac$chr %in% chr_list,]
	ref_h3k27ac <- makeGRangesFromDataFrame(ref_h3k27ac, keep.extra.columns=TRUE)
	

	
	query_h3k27ac <- fread(paste0("C:/Users/pubud/Documents/R_analysis/chip/",mark,"/make peaks/",df))
	

colnames(query_h3k27ac) <- c("chr", "start", "end" )

query_h3k27ac <- query_h3k27ac[query_h3k27ac$chr %in% chr_list,]
	
	query_h3k27ac <- makeGRangesFromDataFrame(query_h3k27ac, keep.extra.columns=TRUE)
		tot <- length(query_h3k27ac)
	overlap_my_ref <- countOverlaps(  query_h3k27ac,ref_h3k27ac, ignore.strand=TRUE,type="any")

	overlap_my_ref <- subset(overlap_my_ref, overlap_my_ref > 0)
overlap <- 	length(overlap_my_ref)
return(overlap/tot*100)
}
	

## ------------------------------------------------------------------------
chip_overlap_barplot_draw <- function(df=NULL, ylab=paste0("% Peak overklap with \nRoadmap")){
  
  ggplot(df, aes(fill=Tissue, y=V3, x=patient)) + 
    geom_bar(position="dodge", stat="identity", colour="black", size=2)+
   theme_bw()+
   scale_fill_manual(values=c("red", "blue"))+
    scale_y_continuous( limit=c(0,100))+
          #  guides(fill=guide_legend(nrow=2,byrow=TRUE))+
          theme(plot.title = element_text(size=12), legend.position = "bottom", text = element_text(size=29),
          axis.text = element_text( size=21), legend.key.size = unit(1.5, 'lines'), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
						panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.ticks.length=unit(.25, "cm"))+labs(x="Patient",y=ylab)
 
}



## ---- fig.height=6, fig.width=6------------------------------------------

#importantn note:
#These chip peak files were made from original peak files obtained from MACS2. Just select the chromosome coordinates from the files and create new files.


files <- list.files(path = "C:/Users/pubud/Documents/R_analysis/chip/H3k27ac/make peaks/", pattern = ".txt", all.files = T,
           full.names = FALSE, recursive = FALSE,
           ignore.case = T)
 patient <- substr(files, 1, 4)
 Tissue <- substr(files, 6, 6)
 
 df <- data.frame(patient, Tissue=as.character(Tissue), stringsAsFactors=FALSE)
 df2 <- data.frame(patient, Tissue=as.character(Tissue), stringsAsFactors=FALSE)
for(i in 1:length(files)){
 
  
  df[i,3] <- chi_overlap_barplot(files[i], mark="H3k27ac")
  df2[i,3] <- chi_overlap_barplot(files[i], mark="H3k27ac", random_ref=T)
}

 df$Tissue[df$Tissue == "N"] <- "Normal"
 df$Tissue[df$Tissue == "T"] <- "Tumor"
 
  df2$Tissue[df2$Tissue == "N"] <- "Normal"
 df2$Tissue[df2$Tissue == "T"] <- "Tumor"
 
 chip_overlap_barplot_draw(df)
 a <- df %>% subset(Tissue=="Normal")  
   mean(a$V3)
   df
   
   
    chip_overlap_barplot_draw(df2, ylab=paste0("% Peak overklap with \nrandom regions"))
 a <- df2 %>% subset(Tissue=="Normal")  
   mean(a$V3)
   df2

## ---- fig.height=6, fig.width=6------------------------------------------
files <- list.files(path = "C:/Users/pubud/Documents/R_analysis/chip/H3k4me1/make peaks/", pattern = ".txt", all.files = T,
           full.names = FALSE, recursive = FALSE,
           ignore.case = T)
 patient <- substr(files, 1, 4)
 Tissue <- substr(files, 6, 6)
 
 df <- data.frame(patient, Tissue=as.character(Tissue), stringsAsFactors=FALSE)
 df2 <- data.frame(patient, Tissue=as.character(Tissue), stringsAsFactors=FALSE)
for(i in 1:length(files)){
 
  
  df[i,3] <- chi_overlap_barplot(files[i], mark="H3k4me1")
  df2[i,3] <- chi_overlap_barplot(files[i], mark="H3k4me1", random_ref=T)
}

 df$Tissue[df$Tissue == "N"] <- "Normal"
 df$Tissue[df$Tissue == "T"] <- "Tumor"
 
  df2$Tissue[df2$Tissue == "N"] <- "Normal"
 df2$Tissue[df2$Tissue == "T"] <- "Tumor"
 
 
 chip_overlap_barplot_draw(df)
 a <- df %>% subset(Tissue=="Normal")  
   mean(a$V3)
   df
   
    chip_overlap_barplot_draw(df2,ylab=paste0("% Peak overklap with \nrandom regions"))
 a <- df2 %>% subset(Tissue=="Normal")  
   mean(a$V3)
   df2


## ---- fig.height=6, fig.width=6------------------------------------------
files <- list.files(path = "C:/Users/pubud/Documents/R_analysis/chip/H3k4me3/make peaks/", pattern = ".txt", all.files = T,
           full.names = FALSE, recursive = FALSE,
           ignore.case = T)
 patient <- substr(files, 1, 4)
 Tissue <- substr(files, 6, 6)
 
 df <- data.frame(patient, Tissue=as.character(Tissue), stringsAsFactors=FALSE)
 df2 <- data.frame(patient, Tissue=as.character(Tissue), stringsAsFactors=FALSE)
for(i in 1:length(files)){
 
  
  df[i,3] <- chi_overlap_barplot(files[i], mark="H3k4me3")
  df2[i,3] <- chi_overlap_barplot(files[i], mark="H3k4me3", random_ref = T)
}

 df$Tissue[df$Tissue == "N"] <- "Normal"
 df$Tissue[df$Tissue == "T"] <- "Tumor"
 
  df2$Tissue[df2$Tissue == "N"] <- "Normal"
 df2$Tissue[df2$Tissue == "T"] <- "Tumor"
 
 chip_overlap_barplot_draw(df)
 a <- df %>% subset(Tissue=="Normal")  
   mean(a$V3)
df

chip_overlap_barplot_draw(df2, ylab=paste0("% Peak overklap with \nrandom regions"))
 a <- df2 %>% subset(Tissue=="Normal")  
   mean(a$V3)
df2

## ------------------------------------------------------------------------
library(gplots)
library(ggplot2)

## ------------------------------------------------------------------------
heatmap.range <- 800
  srange <- 1
  colors = c(seq(0.0001,0.002,length=800))
  my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 799)
mat.final.promoter <- read.table(
  "170716/promoter.final.R354N.H3K27ac.txt",   sep="\t", header=FALSE)
mat.final.genebody <- read.table(
  "170716/genebody.final.R354N.H3K27ac.txt",   sep="\t", header=FALSE)
mat.final.threeprime <- read.table(
  "170716/threeprime.final.R354N.H3K27ac.txt",   sep="\t", header=FALSE)
mat.final <- cbind(mat.final.promoter[,1:800],mat.final.genebody[,1:800],mat.final.threeprime[,1:800])

gplots::heatmap.2(as.matrix(mat.final),Rowv=F,Colv=F,dendrogram = "none", trace="none", col=my_palette, 
                  breaks=colors, labRow = F, ylab = "Gene", 
                  labCol  = c(" ", rep(" ", ((heatmap.range)-srange)), "TSS",rep(" ", ((heatmap.range*2)-2)),
                              " ") ,margins = c(4, 2), cexRow=0.8 ,
                  xlab="Nucleotide position", 
                  main = "R354N",
                  symm=F,symkey=F,symbreaks=T,scale="none")

###tumor
mat.final.promoter <- read.table(
  "170716/promoter.final.R354T.H3K27ac.txt",   sep="\t", header=FALSE)
mat.final.genebody <- read.table(
  "170716/genebody.final.R354T.H3K27ac.txt",   sep="\t", header=FALSE)
mat.final.threeprime <- read.table(
  "170716/threeprime.final.R354T.H3K27ac.txt",   sep="\t", header=FALSE)
mat.final <- cbind(mat.final.promoter[,1:800],mat.final.genebody[,1:800],mat.final.threeprime[,1:800])

gplots::heatmap.2(as.matrix(mat.final),Rowv=F,Colv=F,dendrogram = "none", trace="none", col=my_palette, 
                  breaks=colors, labRow = F, ylab = "Gene", 
                  labCol  = c(" ", rep(" ", ((heatmap.range)-srange)), "TSS",rep(" ", ((heatmap.range*2)-2)),
                              " ") ,margins = c(4, 2), cexRow=0.8 ,
                  xlab="Nucleotide position", 
                  main = "R354T",
                  symm=F,symkey=F,symbreaks=T,scale="none")

###Using the same code generate heatmaps for other samplesa and other histone marks
#how ever you have create the files using and modifying these scripts.
#all files are not attached

## ------------------------------------------------------------------------

#gene list was obtaining when create the heatmap matrix

chip.gene.list <- read.table(
  "gene.list.txt",   sep="\t", header=FALSE)

#this file was created from ordered.down.regulated.genes.csv
#select the gene and base mean
gene.expression.list <- read.csv(
  "basemean.gene.expression.csv", header=TRUE)
gene.expression.list <- gene.expression.list[,c(1,2)]

#load basemean.gene.expression.csv
gene.expression.level <- merge(chip.gene.list, gene.expression.list, by.x="V1", by.y="X")
gene.expression.level <- data.frame(gene.expression.level[order(gene.expression.level$V2),])
#remove rows with NA
gene.expression.level <- gene.expression.level[complete.cases(gene.expression.level),]
cols <- c("blue", "red")[(gene.expression.level$Normal > 0) + 1] 
barplot(rev(log2(gene.expression.level$Normal+1)), col = rev(cols), horiz=TRUE, border = NA)
barplot(rev(gene.expression.level$Normal), col = rev(cols), horiz=TRUE, border = NA)

ggplot(gene.expression.level, aes(x=V1, y=Normal)) +
  geom_bar(stat='identity') +
  coord_flip()

write.csv(gene.expression.level, file="gene.expression.level.and.chip.seq.enrichment.gene.order.csv")


## ------------------------------------------------------------------------

setwd("C:/Users/pubud/Documents/R_analysis/chip/visualize/H3K27ac")
colors = c(seq(-0.3,0.3,length=800))
heatmap.range = 800
srange =1
flanking.size=4000
erange =(flanking.size*2)
#blue neagitive red positive
#blue down regulated red up regulated
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 799)

ordered.down.regulated.genes <- read.csv(
  "ordered.down.regulated.genes.csv", header=TRUE)

#R354
mat.final <- read.table(
  "finalmatrix.R354N.H3K27ac.diff.cn.txt",   sep="\t", header=FALSE)

#filtered by lg2fc SE < 0.5
#in order to do that select the filtered gene list when find the difference in HPC
#mat.final <- read.table("finalmatrix354.filtered.txt",   sep="\t", header=FALSE)

gplots::heatmap.2(as.matrix(mat.final),Rowv=F,Colv=F,dendrogram = "none", trace="none", col=my_palette, 
                  breaks=colors, labRow = F, ylab = "Gene", 
                  labCol  = c("-4000", rep(" ", ((heatmap.range/2)-srange)), "TSS", rep(" ", ((heatmap.range/2)-srange-1)),
                              "+4000") ,margins = c(4, 2), cexRow=0.8 ,
                  xlab="Nucleotide position",lmat = rbind(c(4,3),c(2,1)), lwid = c(1.5,3), lhei = c(1.5,5), 
                  main = "R354",
                  symm=F,symkey=F,symbreaks=T,scale="none")

#making the heat map for log2foild change of gene expression


R354.gene.order.Rmerge <- read.table(
  "ordered.peak.R354N.H3K27ac.gene.order.Rmerge.txt",  sep="\t", header=TRUE)

#select genes only in Rmerge
#write the code

merged.gene.list.lg2fldchnge <- merge(R354.gene.order.Rmerge, ordered.down.regulated.genes , by.x="gene", by.y="X")
merged.ordered.gene.list.lg2fldchnge <- merged.gene.list.lg2fldchnge[order(merged.gene.list.lg2fldchnge$order),]
#make a color column using the log2foldchange column vector
cols <- c("blue", "red")[(merged.ordered.gene.list.lg2fldchnge$log2FoldChange > 0) + 1]  
#draw the bar plot using color values
barplot(rev(merged.ordered.gene.list.lg2fldchnge$log2FoldChange), col = rev(cols), horiz=TRUE, border = NA)
View(data.frame(cols))


