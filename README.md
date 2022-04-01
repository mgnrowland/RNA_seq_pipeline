# RNAseq_Ashwood_2021
```
To examine how activation of different toll-like receptors (TLR) impacts gene expression in Autism Spectrum Disorder (ASD), 
we cultured peripheral blood monocytes from children with ASD, Pervasive Developmental Disorder Not Otherwise Specified 
(PDDNOS) and typically developing children and treated them with either lipoteichoic acid (LTA) or lipopolysaccharide (LPS)
to activate TLR2 or 4 respectively. Following 24 hours of stimulation, we then performed RNA sequencing to profile mRNA 
responses between non-treated (NT), LTA and LPS treated samples for each diagnosis (TD, ASD or PDDNOS).

Samples deposited at (GEO accession #GSE140702; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140702)

```

#run592 /Data/2c7dycd7ut/UnalignedL1/Project_PAAC_L1_H1401P_Ciernia
#run590 /Data/xfm70hwv8e/UnL6L7/Project_PAAC_H1401P_Ciernia
#run590 /Data/49bxhv50b5/UnL6L7/Project_PAAC_H1401P_Ciernia


#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/

#######################################################################################
#load modules
module load bbmap/37.68 

#######################################################################################
#move to the experiment folder
#make the following folders:
#mkdir fastqc_pretrim
#mkdir fastqc_posttrim
#mkdir trimmed_sequences

#place a copy of truseq.fa.gz into the top level folder
#this .fa contains the adapter for trimming

#######################################################################################

#make poly A tail file to filter out 
printf ">polyA\nAAAAAAAAAAAAA\n>polyT\nTTTTTTTTTTTTT\n" | gzip - >  polyA.fa.gz

#make sample files.txt for each sample
#from within fastq file folder: ls -R *fastq.gz > samples.txt
#all samples: find . -name '*fastq.gz' -print > samples.txt

#######################################################################################
#run trim.sh srcipt:
./trim_L1.sh
./trim_L2.sh

#trims and performs pre and post trim qc for each batch of sequencing independently
#collect all qc together with:
module load multiqc/1.6

multiqc fastqc_pretrim/ --filename PreTrim_multiqc_report.html --ignore-samples Undetermined* --interactive

multiqc fastqc_posttrim/ --filename PostTrim_multiqc_report.html --ignore-samples Undetermined* --interactive

#######################################################################################
#concatenate trimmed reads for each sample together into 1 fastq.gz for alignment
#cat 101123NT*.fastq.gz > 101123NT.fastq.gz
./concatenate.sh 
#######################################################################################

#all samples: find . -name '*fastq.gz' -print > samples.txt
#run ./mm10STARbuild.sh to generate START indexes for 100bp reads for ensembl genes
#run ./START_alignmm10.sh to align to mm10

for sample in `cat samples.txt`
do
R1=${sample}

echo ${R1} "unzipping"

gunzip trimmed_sequences/combinedfastq/${sample}.fastq.gz

echo ${R1} "mapping started"

STAR --runThreadN 12 --genomeDir STAR/GRCh38.p12/star_indices/ --readFilesIn  trimmed_sequences/combinedfastq/${sample}.fastq --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI nM AS MD --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix star_out/${sample} 


echo ${R1} "mapping completed"
done

#Indexed bam files are necessary for many visualization and downstream analysis tools
cd star_out
for bamfile in */starAligned.sortedByCoord.out.bam ; do samtools index ${bamfile}; done

#From this point any further analysis can be applied
   
## code for taking counts.txt file outputs from FeatureCounts and importing into R for EdgeR and/or limmavoom analysis
## loops through each txt file and takes the GeneID and accepted_hits.bam counts columns (counts per gene)
## writes a new text file samplename.out.txt that contains only these two columns
## these files are then read into EdgeR 
## a targets dataframe is made that contains a list of the .out.txt files and the group they belong to
## the targets dataframe is then fed into EdgeR using readDGE function of EdgeR
#9.5.19 AVC
#removed patient 101712 in the ASD group has Fragile X
#added NT vs NT comparisons
#remove PDDNOS condtion from analysis and analyze as a separate group
##May 2021 Megan Rowland

library(ggplot2)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(tidyr)
library(readxl)
#loop for combining GeneID and counts for each txt file
path ="C:\\Users\\Megan\\Documents\\Ciernia Lab\\Bioninformatics\\HumanMacrophageAshwood\\htcounts"

setwd(path)

out.file<-""
group<-""
file.names <- dir(path, pattern ='*counts.txt')

for(i in 1:length(file.names)){
  file <- read.table(file.names[i],skip=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  counts <- data.frame(file$Geneid,file[,7])
  #extract file name without count.txt
  m <- as.data.frame(file.names[i])
  names(m) <-  c("split")
  m$split <- as.character(m$split)
  m <-tidyr::separate(m,split,into= c("name","stuff","stuff2"),sep="\\.")
 
  print(m$name)
  #replace column headings with names
  names(counts) <- c("Geneid",m$name)
  
  #write to output file
  out.file <- counts
  write.table(out.file, file =paste(m$name, ".out.txt", sep=""), sep="\t",quote=F)
}

##The following code is for analysis of RNASeq in Limma Voom package
##########################################################################################

#take files from loop above as new input
input.files <-dir(path, pattern ='*out.txt')
input.files <- as.data.frame(input.files)
names(input.files) <- c("filenames")

#get treatments names
input.files$treatment <-as.character(input.files$filenames)
input.files$treatment <-substr(input.files$treatment,7,nchar(input.files$treatment))
input.files$treatment <- gsub(".out.txt","",input.files$treatment)
input.files$subject.name <-substr(input.files$filenames,1,6)
##########################################################################################
#read in subject info
##########################################################################################
library(readxl)
library(dplyr)
path ="C:\\Users\\Megan\\Documents\\Ciernia Lab\\Bioninformatics\\HumanMacrophageAshwood"
setwd(path)

Info <- read.csv("APP Y3 monocytes updated subject list 4-21-21.csv")

#more info
mInfo <- read_excel("/Users/Megan/Documents/Ciernia Lab/Bioninformatics/HumanMacrophageAshwood/APP\ Y3\ Monocyte\ RNA\ inventory\ 11-2018.xlsx",sheet = "Subject details")

Info1 <- merge(Info,mInfo, by.x="subject.name",by.y="Milo_Subjects.IBC",all=T)
Info1 <- Info1[-c(10)]
names(Info1)[3] <- "Year3_Clinical_Impressions"

Info2 <- merge(input.files,Info1,by.x="subject.name",by.y="subject.name",all=T)

#define combined condition
Info2$Group <- paste(Info2$Group.for.analysis,Info2$treatment, sep=".")


allsamples <- table(Info2$Group)
##########################################################################################
#keep PDDNOS as a separate group
#remove no confirmation
#remove.samples <- Info3 %>% filter(Year3_Clinical_Impressions.y == "no confirmation")

#Remove samples:
#101716 no confirmation
#101738 not ASD
#101712 no FMRP expression = FXS patient 
#101150 Brother has autism 6/30/2011

Info2$subject.name <- as.character(Info2$subject.name)
keep.samples <- Info2 %>% 
  filter(subject.name != "101712") %>%
  filter(subject.name != "101150") %>%
  filter(Year3_Clinical_Impressions != "no confirmation") %>%
  filter(Year3_Clinical_Impressions != "Not ASD" )
  
allsamples <- table(keep.samples$Group)

keep.samples$fastq.gz <- paste(keep.samples$subject.name,keep.samples$treatment,".fastq.gz",sep="")

path ="C:\\Users\\Megan\\Documents\\Ciernia Lab\\Bioninformatics\\HumanMacrophageAshwood\\"
setwd(path)

library(xlsx)
write.csv(keep.samples, "KeptSamples.csv")

#write out kept counts files:
path ="/Users/aciernia/Sync/collaborations/Ashwood/HumanMacrophageLPS/Analysis_4_2021/counts_kept"
setwd(path)
input.files <-dir(path, pattern ='*counts.txt')

keep.samples$filenames2 <- gsub("out", "counts", keep.samples$filenames)

input.files.keep <- input.files[which(input.files %in% keep.samples$filenames2)]
 
for(i in 1:length(input.files.keep)){
file <- read.table(input.files.keep[i],skip=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
#get file name
m <- as.data.frame(input.files.keep[i])
names(m) <-  c("split")
m$split <- as.character(m$split)
m <-tidyr::separate(m,split,into= c("name","stuff","stuff2"),sep="\\.")

print(m$name)

write.table(file, file =paste("C:\\Users\\Megan\\Documents\\Ciernia Lab\\Bioninformatics\\HumanMacrophageAshwood\\htcounts",
m$name, ".counts.txt", sep=""), sep="\t",quote=F)
}


#read in only kept counts files:
path ="C:\\Users\\Megan\\Documents\\Ciernia Lab\\Bioninformatics\\HumanMacrophageAshwood\\htcounts"
setwd(path)
#input.files <-dir(path, pattern ='*out.txt')

#input.files.keep <- input.files[which(input.files %in% keep.samples$filenames)]

input.files.keep <- keep.samples$filenames




##########################################################################################
# Make DEGList
library(edgeR)
library(limma)
library(dplyr)

RG <- readDGE(input.files.keep,group=keep.samples$Group)
colnames(RG$counts) <- gsub(".out","",colnames(RG$counts))

#save
path ="C:\\Users\\Megan\\Documents\\Ciernia Lab\\Bioninformatics\\HumanMacrophageAshwood\\test"
setwd(path)
save(RG,keep.samples, file="RG_RNAseqDEGlist.RData")

# number of unique transcripts
dim(RG)
#58735   130


###########################################################
# Filter for one count per million in at least 1/4 of libraries = 130/4 = 32
keep <- rowSums(cpm(RG)>1)>=32
RGcpm <- RG[keep,]
dim(RGcpm)
# 12535   130

geneid <- rownames(RGcpm) #ensemble IDs for biomart

#save
save(RG,RGcpm,file="DEGlist.RData")
##################################################################################
#read in gene info
##################################################################################

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                host = "asia.ensembl.org")
#head(listAttributes(human),100)
#attributes <- listAttributes(mouse)
#attributes[grep("exon", attributes$name),]


genes <- getBM(attributes = c("ensembl_gene_id","chromosome_name","start_position","end_position","hgnc_symbol","entrezgene_id","description"), 
               filter= "ensembl_gene_id", 
               values = rownames(RGcpm),
               mart = human)


genes$genelength <- abs(genes$end_position - genes$start_position) 

#remove duplicates:keeps only first entry
genes <- genes[!duplicated(genes$ensembl_gene_id),]

#match order: df[match(target, df$name),]
genes <- genes[match(geneid,genes$ensembl_gene_id),]

#save(genes,file="Hg38_Ensemble_geneInfo.RData")
#setwd("/Users/annieciernia/Box/Ashwood_labcollaboration/Experiments/HumanMacrophageLPS/LexogenRNAseq/7_29_19_analysis")
#load("Hg38_Ensemble_geneInfo.RData")

##########################################################################
#add gene info to DEGlist object
##########################################################################

#add into DEGlist
RGcpm$genes <- genes

#save
save(RG,RGcpm,file="DEGlist.RData")

###########################################################
#filtering plot
###########################################################
library(RColorBrewer)
nsamples <- ncol(RGcpm)

colourCount = nsamples
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fill=getPalette(colourCount)

#plot:
pdf('FilteringCPM_plots.pdf')
par(mfrow=c(1,2))

#prefilter:
lcpm <- cpm(RG, log=TRUE, prior.count=2)


plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")

#filtered data
#og-CPM of zero threshold (equivalent to a CPM value of 1) used in the filtering ste
lcpm <- cpm(RGcpm, log=TRUE, prior.count=2)
plot(density(lcpm[,1]), col=fill[1], lwd=2, ylim=c(0,0.5), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=fill[i], lwd=2)
}
#legend("topright", Samples, text.col=fill, bty="n")
dev.off()


###########################################################
###########################################################
#reset library sizes
RGcpm$samples$lib.size <- colSums(RGcpm$counts)


#plot library sizes
pdf('LibrarySizes.pdf',w=30,h=8)
barplot(RGcpm$samples$lib.size,names=colnames(RGcpm),las=2)
# Add a title to the plot
title("Barplot of library sizes")
dev.off()


# Get log2 counts per million
logcounts <- cpm(RGcpm,log=TRUE)
# Check distributions of samples using boxplots
pdf('NonNormalizedLogCPM.pdf',w=30,h=10)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()

##########################################################################
#Diagnosis * treatment + sex with subjects correlation removed
##################################################################################
#make design matrix
keep.samples$sex <- factor(keep.samples$`SUBJ GENDER`, levels=c("1","2"))

keep.samples$treatment <- factor(keep.samples$treatment,levels = c("NT","LPS","LTA"))

keep.samples$Diagnosis <- keep.samples$Group.for.analysis

keep.samples$Diagnosis <- factor(keep.samples$Diagnosis, levels = c("TD","ASD","PDDNOS"))

keep.samples$Group <- paste(keep.samples$Diagnosis,keep.samples$treatment,sep=".")
keep.samples$Group <- factor(keep.samples$Group)

#set design matrix
design <- model.matrix(~0+Group+sex, data=keep.samples) #if using a 0 intercept must set up contrasts

colnames(design)

#make contrasts
#https://support.bioconductor.org/p/91718/
#test whether the average treatment effect across all cell lines and treatment methods is significantly different from zero
#contrast compares the average of all treatment samples to the average of all control samples.

#https://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
#Contrasts NT vs LPS and NT vs LTA for each diagnosis:
my.contrasts <- makeContrasts(
  #NT vs NT 
  Typical.NTvsASD.NT = GroupASD.NT - GroupTD.NT,
  Typical.NTvsPDDNOS.NT =  GroupPDDNOS.NT - GroupTD.NT,
  Typical.NTvsTD.LPS =  GroupTD.LPS - GroupTD.NT,
  Typical.NTvsTD.LTA =  GroupTD.LTA - GroupTD.NT,
  ASD.NTvsASD.LPS =  GroupASD.LPS - GroupASD.NT,
  ASD.NTvsASD.LTA =  GroupASD.LTA - GroupASD.NT,
  PDDNOS.NTvsASD.NT = GroupASD.NT - GroupPDDNOS.NT,
  PDDNOS.NTvsPDDNOS.LPS =  GroupPDDNOS.LPS - GroupPDDNOS.NT,
  PDDNOS.NTvsPDDNOS.LTA =  GroupPDDNOS.LTA - GroupPDDNOS.NT,
  levels=design)
  #interactions: NT vs LPS for ASD vs Typical:
  #ASDNTvsLPS_TypicalNTvsLPS = (GroupASD.LPS - GroupASD.NT) - (GroupTypical.LPS-GroupTypical.NT),
  #interactions: NT vs LPS for PDDNOS vs Typical:
  #ASDNTvsLPS_TypicalNTvsLPS = (GroupPDDNOS.LPS - GroupPDDNOS.NT) - (GroupTypical.LPS-GroupTypical.NT),
  #interactions: NT vs LPS for PDDNOS vs ASD:
  #ASDNTvsLPS_TypicalNTvsLPS = (GroupPDDNOS.LPS - GroupPDDNOS.NT) - (GroupASD.LPS-GroupASD.NT),
  #interactions: NT vs LTA for ASD vs Typical:
 # ASDNTvsLTA_TypicalNTvsLTA = (GroupASD.LTA - GroupASD.NT) - (GroupTypical.LTA-GroupTypical.NT),
  #interactions: NT vs LTA for PDDNOS vs Typical:
 # ASDNTvsLTA_TypicalNTvsLTA = (GroupPDDNOS.LTA - GroupPDDNOS.NT) - (GroupTypical.LTA-GroupTypical.NT),
  #interactions: NT vs LTA for PDDNOS vs ASD:
 # ASDNTvsLTA_TypicalNTvsLTA = (GroupPDDNOS.LTA - GroupPDDNOS.NT) - (GroupASD.LTA-GroupASD.NT),


##################################################################################
#Normalization TMM and Voom
##################################################################################
#The two steps refer to different aspects of normalization. CPM "normalization" accounts for library size differences between samples, and produces normalized values that can be compared on an absolute scale (e.g., for filtering). TMM normalization accounts for composition bias, and computes normalization factors for comparing between libraries on a relative scale. CPM normalization doesn't account for composition bias, and TMM normalization doesn't produce normalized values. Thus, you need both steps in the analysis pipeline. This isn't a problem, as the two steps aren't really redundant.
#TMM normalization for library composition
DGE=calcNormFactors(RGcpm,method =c("TMM")) 

pdf('VoomTrend.pdf',w=6,h=4)
v=voom(DGE,design,plot=T)
dev.off()

corfit <- duplicateCorrelation(v, design, block=keep.samples$subject.name)

fit <- lmFit(v, design, block = keep.samples$subject.name, correlation = corfit$consensus)

fit2 <- contrasts.fit(fit,my.contrasts)
fit2 <- eBayes(fit2,robust=TRUE )

pdf('PlotSA_VoomTrend.pdf',w=6,h=4)
plotSA(fit2, main="Final model: Mean variance trend")
dev.off()

#save
save(RG,RGcpm,keep.samples,design,my.contrasts,fit,fit2,v,file="DEGlist.RData")




#box plots for the voom normalised data to compare to before normalisation (only RLE)
#v$E already log2CPM
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="TMM & Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")

#plot only normalized data
pdf('TMM_VOOM_NormalizedLogCPM.pdf',w=30,h=10)
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="TMM & Voom transformed logCPM")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
dev.off()


##########################################################################
#MAA plots
##########################################################################

#MAA plot for each sample
keep.samples$Name <- paste(keep.samples$subject.name,keep.samples$Group,sep=".")
keep.samples$Name <- gsub("\\."," ", keep.samples$Name)

C <- unique(keep.samples$Name)
#Mean-Difference Plot of Expression Data
#8 per page
pdf('Samples_MAAplots.pdf')
par(mfrow=c(4,2),mar=c(5,6,4,2)+0.1)
for (i in 1:length(C)) {
  print(paste(i))
  name <- C[i]
  #pdf(file = paste(name,"MAAplot.pdf", sep ="_"), wi = 3, he = 3)
  #MAA: should center on zero
  plotMD(v,column=i,
         ylab = "Expression log-ratio\n(this sample vs others)")
  abline(h=0, col="red", lty=2, lwd=2)
  #dev.off()
  
}

dev.off()

##################################################################################

pdf(file = "MDSplot_Diagnosis.pdf", wi = 12, he = 10, useDingbats=F)

levels(keep.samples$Diagnosis)
col.cell <- c("blue","orange","green")[keep.samples$Diagnosis]
data.frame(keep.samples$Diagnosis,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(v,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(4, 0.5,
       fill=c("blue","orange","green"),
       legend=levels(keep.samples$Diagnosis),
       cex = 0.8)
# Add a title
title("Diagnosis by Treatment MDS Plot")
dev.off()

#plot MDS for sex
pdf(file = "MDSplot_sex.pdf", wi = 12, he = 10, useDingbats=F)

levels(keep.samples$sex)
col.cell2 <- c("blue","orange")[keep.samples$sex]
data.frame(keep.samples$sex,col.cell2)


# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(v,col=col.cell2)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(4, 0.5,
       fill=c("blue","orange"),
       legend=levels(keep.samples$sex),
       cex = 0.8)
# Add a title
title("Sex MDS Plot")

dev.off()


#by Treatment:
pdf(file = "MDSplot_Treatment.pdf", wi = 12, he = 10, useDingbats=F)

#par(mfrow=c(1,2))
#plot MDS for treatment
levels(keep.samples$treatment)
col.cell <- c("blue","orange","green")[keep.samples$treatment]
data.frame(keep.samples$treatment,col.cell)

# MDS with group colouring
par(xpd = T, mar = par()$mar + c(0,0,0,7))
plotMDS(v,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend(4, 0.5,
       fill=c("blue","orange","green"),
       legend=levels(keep.samples$treatment),
       cex = 0.8)
# Add a title
title("Treatment MDS Plot")
dev.off()

##################################################################################
#Average CPM for each condition
##################################################################################

#get log2CPM counts from voom and put in dataframe:
library(plotrix)
library(xlsx)
#average log2 CPM and sem
countdf <- as.data.frame(v$E)
countdf$GeneID <- rownames(v$E)
#targets <- v$targets
#targets$sample <- gsub(".out","",rownames(targets))

#DF <- merge(countdf,comp_out, by.x ="GeneID",by.y="genes")

#summarize 
countdf2 <- countdf %>% group_by(GeneID) %>% gather(sample,log2CPM, 1:(ncol(countdf)-1)) 
countdf2 <- as.data.frame(countdf2)
keep.samples$sample <- gsub(".out.txt","",keep.samples$filenames)
countdf3 <-merge(countdf2,keep.samples,by="sample")

GeneSummary <- countdf3 %>% group_by(GeneID,Group) %>% 
  summarize(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
 # dplyr::select(GeneID,Group,meanlog2CPM) %>%
  spread(Group,meanlog2CPM)

#add in gene symbol information
geneinfo <- v$genes
GeneSummary2 <- merge(GeneSummary,geneinfo ,by.x="GeneID",by.y="ensembl_gene_id")

write.csv(GeneSummary2,file= "MeanLog2CPM_pergroup_GeneExpression.csv")

#log2CPM values for individuals
countdf3$Individual <- paste(countdf3$Group,countdf3$sample, sep=".")
countdf3$Individual <- factor(countdf3$Individual)

IndGeneSummary <- countdf3 %>% dplyr::group_by(GeneID,Individual) %>% 
  summarise(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
  # dplyr::select(GeneID,Group,meanlog2CPM) %>%
  spread(Individual,meanlog2CPM)

#add in gene symbol information
geneinfo <- v$genes
IndGeneSummary2 <- merge(IndGeneSummary,geneinfo ,by.x="GeneID",by.y="ensembl_gene_id")

write.csv(IndGeneSummary2,file= "IndividualSamplesLog2CPM_GeneExpression.csv")


#go back to CPM from log2: 2^log2CPM = CPM
GeneSummaryCPM <- countdf3 %>% group_by(GeneID,Group) %>% 
  summarize(meanlog2CPM = mean(log2CPM)) %>%
  ungroup()  %>%
  mutate(meanCPM = 2^meanlog2CPM) %>%
  dplyr::select(GeneID,Group,meanCPM) %>%
  spread(Group,meanCPM)

#add in gene symbol information
GeneSummaryCPM2 <- merge(GeneSummaryCPM,geneinfo,by.x="GeneID",by.y="ensembl_gene_id")

write.csv(GeneSummaryCPM2,file= "CPM_GeneExpression.csv")

save(GeneSummary2,GeneSummaryCPM2,file="AverageCPMandLog2CPM.RData")

##################################################################################
#DE analysis:
##################################################################################
str(fit2)
#summary(decideTests(fit2,adjust.method="fdr", method="global"))
#sum <- summary(decideTests(fit2,adjust.method="BH", method="separate"))

#more weight to fold-changes in the ranking
#treat computes empirical Bayes moderated-t p-values relative to a minimum meaningful fold-change threshold. 
#Instead of testing for genes that have true log-fold-changes different from zero, it tests whether the true 
#log2-fold-change is greater than lfc in absolute value (McCarthy and Smyth, 2009). 
#In other words, it uses an interval null hypothesis, where the interval is [-lfc,lfc]. 
#When the number of DE genes is large, treat is often useful for giving preference to larger fold-changes and 
#for prioritizing genes that are biologically important. 
#treat is concerned with p-values rather than posterior odds, so it does not compute the B-statistic lods. 
#The idea of thresholding doesn't apply to F-statistics in a straightforward way, so moderated F-statistics are also not computed. 
#When lfc=0, treat is identical to eBayes, except that F-statistics and B-statistics are not computed. 
#The lfc threshold is usually chosen relatively small, because significantly DE genes must all have fold changes substantially greater than the testing threshold. 
#Typical values for lfc are log2(1.1), log2(2) or log2(2). The top genes chosen by treat can be examined using topTreat.

# Treat relative to a ~2 fold-change
tfit <- treat(fit2,lfc=log2(2))
dt <- decideTests(tfit,adjust.method="BH", method="separate")
sum <- summary(dt)
sum
write.csv(sum,"SummaryCount_DEGs.csv")

dt <- as.data.frame(dt)


#######################################################################
#get out DE lists for each contrast:
#######################################################################

comparisons=(coef(tfit))
comparisons=colnames(comparisons)

comp_out <- as.data.frame(rownames(v$E))
names(comp_out) <- c("GeneID")
nrowkeep <- nrow(comp_out)

SumTableOut <- NULL


for(i in 1:length(comparisons)){
  #comparison name
  comp=comparisons[i]
  print(comp)
  #make comparisons 
  
  tmp=topTreat(tfit,coef=i,number=nrowkeep,adjust.method="BH")
  #nrow(tmp[(tmp$adj.P.Val<0.05),]) # number of DE genes
  
  #LogFC values:https://support.bioconductor.org/p/82478/
  tmp$direction <- c("none")
  tmp$direction[which(tmp$logFC > 0)] = c("Increase")
  tmp$direction[which(tmp$logFC < 0)] = c("Decrease")
  
  tmp$significance <- c("nonDE")
  tmp$significance[which(tmp$adj.P.Val <0.05)] = c("DE")
  
  #summary counts table based on Ensemble Gene ID counts:
  SumTable <- table(tmp$significance,tmp$direction)
  SumTable <- as.data.frame(SumTable)
  SumTable$comparison <- paste(comp)
  SumTableOut <- rbind(SumTable,SumTableOut)
  
  #get geneids  
  tmp$GeneID <- rownames(tmp)
  
  #gene gene names and expression levels
  tmp2 <- tmp
  
  tmp2$comparison <- paste(comp)
  
  write.csv(tmp2,file = paste(comp,"_DEgenes.csv"))
  
  #save to output:
  #merge <- merge(comp_out,tmp2, by= "GeneID")
  merge2 <- tmp2 %>% dplyr::select(ensembl_gene_id,logFC,t,P.Value,adj.P.Val,direction,significance)
  colnames(merge2) <- paste(colnames(merge2),comp,sep=".")
  colnames(merge2)[1] <- c("GeneID")
  comp_out <- merge(comp_out, merge2, by="GeneID")

  
  #data for plot with gene names:
  genenames <- tmp2 %>% dplyr::select(adj.P.Val,logFC,hgnc_symbol) %>% distinct()
  
  #names for plots
  plotname <- gsub("\\."," ",comp)
  plotname <- gsub("vs"," vs ",plotname)
  
  #volcano plot
  pdf(file = paste(comp,"_Volcano.pdf", sep=""), wi = 9, he = 6, useDingbats=F)
  
  with(genenames, plot(logFC, -log10(adj.P.Val), pch=20,col="gray", main=paste(plotname,"\nVolcano plot", sep=" "), ylab =c("-log10(adj.pvalue)"),xlab =c("Log Fold Change") ))
  
  #color points red when sig and log2 FC > 2 and blue if log2 FC < -2 
  with(subset(genenames, logFC < -2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  with(subset(genenames, logFC > 2 & -log10(adj.P.Val) > -log10(.05)), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  #add lines
  abline(h = -log10(.05), col = c("black"), lty = 2, lwd = 1)
  abline(v = c(-2,2), col = "black", lty = 2, lwd = 1)
  
  #Label points with the textxy function from the calibrate plot
  library(calibrate)
  with(subset(genenames, adj.P.Val<0.05 & abs(logFC)>10), textxy(logFC, -log10(adj.P.Val), labs=hgnc_symbol, cex=.4))
  
  dev.off()
  
}


write.csv(SumTableOut,"SummaryTableDEgenes.csv")

#master outfile to get log2CPM values
mout <- merge(comp_out, GeneSummary2, by="GeneID")
write.csv(mout,"AllDEG_AllConditions_log2CPM.csv")

save(v,corfit,fit2,tfit,dt,SumTableOut,comp_out,GeneSummaryCPM2,GeneSummary2,geneinfo,countdf,RG,genes,file = "AshwoodRNAseq.RData")
#load("AshwoodRNAseq.RData")
#######################################################################
#DE plots for Upregulated genes with LPS treatment
#######################################################################

#load("AshwoodRNAseq.RData")
#venn for LPS:
#The number of genes that are not DE in either comparison are marked in the bottom-right.

pdf("Venn_LPSgreaterthanNT.pdf",height=10,width = 10)
vennDiagram(dt[,c("Typical.NTvsTD.LPS","ASD.NTvsASD.LPS","PDDNOS.NTvsPDDNOS.LPS")],
            circle.col=c("salmon","purple","grey"), include="down",
            names=c("Typical NT>LPS","ASD NT>LPS","PDDNOS NT>LPS"))
dev.off()

pdf("Venn_LPSlessthanNT.pdf",height=10,width = 10)
vennDiagram(dt[,c("Typical.NTvsTD.LPS","ASD.NTvsASD.LPS","PDDNOS.NTvsPDDNOS.LPS")],
            circle.col=c("salmon","purple","grey"), include="up",
            names=c("Typical NT<LPS","ASD NT<LPS","PDDNOS NT<LPS"))
dev.off()


#######euler
library(eulerr)
fit <- euler(c("A" = 215, "B" = 134, "C" = 48, "A&B" = 260, "A&C" = 34, "B&C" = 31, "A&B&C" = 754))
pdf("euler_LPSgreaterthanNT.pdf",height=10,width = 10)
plot(fit, quantities = TRUE, edges = FALSE,  fills = c("cadetblue4", "brown3", "darkgoldenrod1"), labels = c("Typical NT>LPS","ASD NT>LPS","PDDNOS NT>LPS")) 
dev.off()

fit <- euler(c("A" = 122, "B" = 78, "C" = 34, "A&B" = 248, "A&C" = 25, "B&C" = 13, "A&B&C" = 956))
pdf("euler_LPSlessthanNT.pdf",height=10,width = 10)
plot(fit, quantities = TRUE, edges = FALSE,  fills = c("cadetblue4", "brown3", "darkgoldenrod1"), labels = c("Typical NT<LPS","ASD NT<LPS","PDDNOS NT<LPS")) 
dev.off()


#list of conditions to overlap
Typical.NTlessLPS <- mout %>% filter(significance.Typical.NTvsTD.LPS == "DE" & direction.Typical.NTvsTD.LPS == "Increase")
Typical.NTgreaterLPS <- mout %>% filter(significance.Typical.NTvsTD.LPS == "DE" & direction.Typical.NTvsTD.LPS == "Decrease")
ASD.NTlessLPS <- mout %>% filter(significance.ASD.NTvsASD.LPS == "DE" & direction.ASD.NTvsASD.LPS == "Increase")
ASD.NTgreaterLPS <- mout %>% filter(significance.ASD.NTvsASD.LPS == "DE" & direction.ASD.NTvsASD.LPS == "Decrease")
PDDNOS.NTlessLPS <- mout %>% filter(significance.PDDNOS.NTvsPDDNOS.LPS == "DE" & direction.PDDNOS.NTvsPDDNOS.LPS == "Increase")
PDDNOS.NTgreaterLPS <- mout %>% filter(significance.PDDNOS.NTvsPDDNOS.LPS == "DE" & direction.PDDNOS.NTvsPDDNOS.LPS == "Decrease")

#combine into list
x = list(Typical.NTgreaterLPS= Typical.NTgreaterLPS$GeneID,
         ASD.NTgreaterLPS= ASD.NTgreaterLPS$GeneID,
         PDDNOS.NTgreaterLPS= PDDNOS.NTgreaterLPS$GeneID,
         Typical.NTlessLPS=Typical.NTlessLPS$GeneID,
         ASD.NTlessLPS=  ASD.NTlessLPS$GeneID,
         PDDNOS.NTlessLPS=  PDDNOS.NTlessLPS$GeneID)


rapply(x, length, how="list")

#$Typical.NTgreaterLPS
#[1] 1261

#$ASD.NTgreaterLPS
#[1] 1177

#$PDDNOS.NTgreaterLPS
#[1] 865

#$Typical.NTlessLPS
#[1] 1349

#$ASD.NTlessLPS
#[1] 1293

#$PDDNOS.NTlessLPS
# 1027

##################################################################
#Fisher's exact test for list overlaps LPS
##################################################################
library(GeneOverlap)
#https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf

#background = all genes expressed in the experiment
BG <- length(GeneSummaryCPM$GeneID)

#overlap list all vs all
go.obj <-  newGOM(x,genome.size=BG)
print(go.obj)
drawHeatmap(go.obj)

#get p values
options(digits=10)
pval <- getMatrix(go.obj, name="pval")

pval <- as.data.frame(pval)
pval_df <- pval %>% mutate(List1 = rownames(pval)) %>%
  gather(List2,FisherExactpvalue,1:3)

#get odds ratio
OR <- getMatrix(go.obj, "odds.ratio")
OR <- as.data.frame(OR)

OR_df <- OR %>% mutate(List1 = rownames(OR)) %>%
  gather(List2,OddsRatio,1:3)

Overlaps <- merge(pval_df,OR_df,by=c("List1","List2"))

#remove List1 == List2
Overlaps <- Overlaps[Overlaps$List1!=Overlaps$List2,]

#adjust by fdr
Overlaps$fdr <- p.adjust(as.numeric(Overlaps$FisherExactpvalue),method="fdr")

write.csv(Overlaps,"LPSVennListOverlapsFishersExactTest.csv")

##################################################################
#plot:
##################################################################

library(nVennR)

#get overlapping genes:
myV <- plotVenn(x,showPlot = T)
listregions <- listVennRegions(myV)

#fix names:
namesregionlist <- names(listregions)

library(qdapRegex)
names(listregions) <- rm_between(namesregionlist, "(", ")", extract=TRUE)

library(tidyverse)
l_tib <- listregions %>% 
  #unlist(recursive = T) %>% 
  enframe() %>% 
  unnest()

l_tib <- as.data.frame(l_tib)
names(l_tib) <- c("LPSVennGroup","GeneID")
l_tib$GeneID <- as.character(l_tib$GeneID)

mout2 <- merge(mout,l_tib, by = "GeneID",all.x=T)

write.csv(mout2,"AllDEG_AllConditions.csv")


#######################################################################
#DE plots for Upregulated genes with LTA treatment
#######################################################################

pdf(file = "Venn_LTAgreaterthanNT.pdf", wi=10, he=10)
vennDiagram(dt[,c("Typical.NTvsTD.LTA","ASD.NTvsASD.LTA", "PDDNOS.NTvsPDDNOS.LTA")],
            circle.col=c("salmon","purple","grey"), include="down",
            names=c("Typical NT>LTA","ASD NT>LTA", "PDDNOS NT>LTA"))

dev.off()

pdf("Venn_LTAlessthanNT.pdf",height=10,width = 10)
vennDiagram(dt[,c("Typical.NTvsTD.LTA","ASD.NTvsASD.LTA","PDDNOS.NTvsPDDNOS.LTA")],
            circle.col=c("salmon","purple","grey"), include="up",
            names=c("Typical NT<LTA","ASD NT<LTA","PDDNOS NT<LTA"))
dev.off()

#######euler

fit <- euler(c("A" = 278, "B" = 94, "C" = 14, "A&B" = 423, "A&C" = 23, "B&C" = 7, "A&B&C" = 369))
pdf("euler_LTAgreaterthanNT.pdf",height=10,width = 10)
plot(fit, quantities = TRUE, edges = FALSE,  fills = c("cadetblue4", "brown3", "darkgoldenrod1"), labels = c("Typical NT>LTA","ASD NT>LTA","PDDNOS NT>LTA")) 
dev.off()

fit <- euler(c("A" = 128, "B" = 86, "C" = 21, "A&B" = 376, "A&C" = 7, "B&C" = 6, "A&B&C" = 645))
pdf("euler_LTAlessthanNT.pdf",height=10,width = 10)
plot(fit, quantities = TRUE, edges = FALSE,  fills = c("cadetblue4", "brown3", "darkgoldenrod1"), labels = c("Typical NT<LTA","ASD NT<LTA","PDDNOS NT<LTA")) 
dev.off()


#list of conditions to overlap
Typical.NTlessLTA <- mout %>% filter(significance.Typical.NTvsTD.LTA == "DE" & direction.Typical.NTvsTD.LTA == "Increase")
Typical.NTgreaterLTA <- mout %>% filter(significance.Typical.NTvsTD.LTA == "DE" & direction.Typical.NTvsTD.LTA == "Decrease")
ASD.NTlessLTA <- mout %>% filter(significance.ASD.NTvsASD.LTA == "DE" & direction.ASD.NTvsASD.LTA == "Increase")
ASD.NTgreaterLTA <- mout %>% filter(significance.ASD.NTvsASD.LTA == "DE" & direction.ASD.NTvsASD.LTA == "Decrease")
PDDNOS.NTlessLTA <- mout %>% filter(significance.PDDNOS.NTvsPDDNOS.LTA == "DE" & direction.PDDNOS.NTvsPDDNOS.LTA == "Increase")
PDDNOS.NTgreaterLTA <- mout %>% filter(significance.PDDNOS.NTvsPDDNOS.LTA == "DE" & direction.PDDNOS.NTvsPDDNOS.LTA == "Decrease")

#combine into list
x2 = list(Typical.NTgreaterLTA= Typical.NTgreaterLTA$GeneID,
         ASD.NTgreaterLTA= ASD.NTgreaterLTA$GeneID,
         PDDNOS.NTlessLTA=  PDDNOS.NTlessLTA$GeneID,
         Typical.NTlessLTA=Typical.NTlessLTA$GeneID,
         ASD.NTlessLTA=  ASD.NTlessLTA$GeneID,
         PDDNOS.NTgreaterLTA=  PDDNOS.NTgreaterLTA$GeneID)


rapply(x2, length, how="list")

#$Typical.NTgreaterLTA
#[1] 1091

#$ASD.NTgreaterLTA
#[1] 890

#$PDDNOS.NTlessLTA
#[1] 679

#$Typical.NTlessLTA
#[1] 1154

#$ASD.NTlessLTA
#[1] 1111

#$PDDNOS.NTgreaterLTA
#[1] 413

##################################################################
#Fisher's exact test for list overlaps LTA
##################################################################
library(GeneOverlap)
#https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf

#background = all genes expressed in the experiment
BG <- length(GeneSummaryCPM$GeneID)

#overlap list all vs all
go.obj <-  newGOM(x2,genome.size=BG)
print(go.obj)
drawHeatmap(go.obj)

#get p values
options(digits=10)
pval <- getMatrix(go.obj, name="pval")

pval <- as.data.frame(pval)
pval_df <- pval %>% mutate(List1 = rownames(pval)) %>%
  gather(List2,FisherExactpvalue,1:3)

#get odds ratio
OR <- getMatrix(go.obj, "odds.ratio")
OR <- as.data.frame(OR)

OR_df <- OR %>% mutate(List1 = rownames(OR)) %>%
  gather(List2,OddsRatio,1:3)

Overlaps <- merge(pval_df,OR_df,by=c("List1","List2"))

#remove List1 == List2
Overlaps <- Overlaps[Overlaps$List1!=Overlaps$List2,]

#adjust by fdr
Overlaps$fdr <- p.adjust(as.numeric(Overlaps$FisherExactpvalue),method="fdr")

write.csv(Overlaps,"LTAVennListOverlapsFishersExactTest.csv")

##################################################################
#plot:
##################################################################


library(nVennR)

#get overlapping genes:
myV <- plotVenn(x2,showPlot = F)
listregions <- listVennRegions(myV)

#fix names:
namesregionlist <- names(listregions)

library(qdapRegex)
names(listregions) <- rm_between(namesregionlist, "(", ")", extract=TRUE)

library(tidyverse)
l_tib <- listregions %>% 
  #unlist(recursive = T) %>% 
  enframe() %>% 
  unnest()

l_tib <- as.data.frame(l_tib)
names(l_tib) <- c("LTAVennGroup","GeneID")
l_tib$GeneID <- as.character(l_tib$GeneID)

mout3 <- merge(mout2,l_tib, by = "GeneID",all.x=T)

write.csv(mout3,"AllDEG_AllConditions.csv")

#merge into individual files:
IndGeneSummary2 <- merge(mout3, IndGeneSummary2, by="GeneID")
write.csv(IndGeneSummary2, file="Individual_log2CPM_allinfo.csv")
##################################################################
#LPS vs LTA: all list comparisons Upset plots
##################################################################

library("UpSetR")

#list of conditions to overlap
listInput = list("Typical.NTvsTypical.LPS" = rownames(dt[which(dt[,"Typical.NTvsTD.LPS"] ==1),]),
                 "Typical.NTvsTypical.LTA" = rownames(dt[which(dt[,"Typical.NTvsTD.LTA"] ==1),]),
                 "ASD.NTvsASD.LPS" = rownames(dt[which(dt[,"ASD.NTvsASD.LPS"] ==1),]),
                 "ASD.NTvsASD.LTA" = rownames(dt[which(dt[,"ASD.NTvsASD.LTA"] ==1),]),
                 "PDDNOS.NTvsPDDNOS.LPS" = rownames(dt[which(dt[,"PDDNOS.NTvsPDDNOS.LPS"] ==1),]),
                 "PDDNOS.NTvsPDDNOS.LTA" = rownames(dt[which(dt[,"PDDNOS.NTvsPDDNOS.LTA"] ==1),])
)

rapply(listInput, length, how="list")

#colors:
queries = 
  list( #only in NT
       list(query =  intersects, params = list("Typical.NTvsTypical.LPS"), color = "orange", active = T,query.name = c("NT only")), 
       list(query =  intersects, params = list("Typical.NTvsTypical.LTA"), color = "orange", active = T,query.name = c("NT only")), 
       list(query =  intersects, params = list("Typical.NTvsTypical.LPS", "Typical.NTvsTypical.LTA"), color = "orange", active = T,query.name = c("NT only")), 

       #only in ASD LPS or LTA
      list(query = intersects,  params = list("ASD.NTvsASD.LPS","ASD.NTvsASD.LTA"), color = "red", active = T,query.name = c("ASD only")), 
      list(query = intersects,  params = list("ASD.NTvsASD.LPS"), color = "red", active = T), #only in ASD LPS or LTA
      list(query = intersects,  params = list("ASD.NTvsASD.LTA"), color = "red", active = T), #only in ASD LPS or LTA

      #only in PDDNOS LPS or LTA
      list(query = intersects,  params = list("PDDNOS.NTvsPDDNOS.LPS","PDDNOS.NTvsPDDNOS.LTA"), color = "blue", active = T,query.name = c("PDDNOS only")), 
      list(query = intersects,  params = list("PDDNOS.NTvsPDDNOS.LPS"), color = "blue", active = T), #only in PDDNOS LPS or LTA
      list(query = intersects,  params = list("PDDNOS.NTvsPDDNOS.LTA"), color = "blue", active = T), #only in PDDNOS LPS or LTA
      
      #only in ASD or PDDNOS LPS or LTA
      list(query = intersects,  params = list("PDDNOS.NTvsPDDNOS.LPS","ASD.NTvsASD.LPS"), color = "pink", active = T,query.name = c("ASD or PDDNOS only")), 
      list(query = intersects,  params = list("PDDNOS.NTvsPDDNOS.LTA","ASD.NTvsASD.LTA"), color = "pink", active = T) #only in ASDorPDDNOS LPS
     
)

#blue = genes regulated only in PDDNOS
#red = genes regulated only in ASD
#orange = genes regulated only in NT
#magenta = genes only in ASD or PDDNOS

pdf(file = "UpSetPlot_AllDEGs.pdf", wi = 14, he = 10, useDingbats=F)
upset(fromList(listInput),sets = rev(names(listInput)), 
      order.by = "freq",
      #number.angles =10,
      keep.order = TRUE,
      text.scale = c(2),
      query.legend = "top",
      queries = queries)

dev.off()

pdf(file = "UpSetPlot_AllDEGs_nocolour.pdf", wi = 14, he = 10, useDingbats=F)
upset(fromList(listInput),sets = rev(names(listInput)), 
      order.by = "freq",
      #number.angles =10,
      keep.order = TRUE,
      text.scale = c(2),
      query.legend = "top"
      #queries = queries
      )

dev.off()

table(mout3$LPSVennGroup,mout3$LTAVennGroup)

##################################################################
#Fisher's exact test for list overlaps LPS
##################################################################
library(GeneOverlap)
#https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf

#background = all genes expressed in the experiment
BG <- length(GeneSummaryCPM$GeneID)

go.obj <-  newGOM(listInput,genome.size=BG)
print(go.obj)

drawHeatmap(go.obj,adj.p = T)

#get p values
options(digits=10)
pval <- getMatrix(go.obj, name=c("pval") )

pval <- as.data.frame(pval)
pval$FisherExactTest <- c("pvalues")

#shared gene counts
intersection <- getMatrix(go.obj,name=c("intersection")) 
intersection <- as.data.frame(intersection)
intersection$FisherExactTest <- c("shared genes")


#get odds ratio
OR <- getMatrix(go.obj, "odds.ratio")
OR <- as.data.frame(OR)
OR$FisherExactTest <- c("odds ratio")

Overlaps <- rbind(pval,OR,intersection)

write.csv(Overlaps,"AllDEGs_VennListOverlapsFishersExactTest.csv")



#######################################################################
# GO enrichment of each list using ClusterCompare
#######################################################################
#library(devtools)
#devtools::install_github("YuLab-SMU/clusterProfiler.dplyr")
library("clusterProfiler")
#The input for geneCluster parameter should be a named list of gene IDs. 
rownames(mout3) <- mout3$GeneID
listDEGs = list("Typical.NTvsTypical.LPS" = rownames(mout3[which(mout3[,"significance.Typical.NTvsTD.LPS"] =="DE"),]),
                 "Typical.NTvsTypical.LTA" = rownames(mout3[which(mout3[,"significance.Typical.NTvsTD.LTA"] =="DE"),]),
                 "ASD.NTvsASD.LPS" = rownames(mout3[which(mout3[,"significance.ASD.NTvsASD.LPS"] =="DE"),]),
                 "ASD.NTvsASD.LTA" = rownames(mout3[which(mout3[,"significance.ASD.NTvsASD.LTA"] =="DE"),]),
                 "PDDNOS.NTvsPDDNOS.LPS" = rownames(mout3[which(mout3[,"significance.PDDNOS.NTvsPDDNOS.LPS"] =="DE"),]),
                 "PDDNOS.NTvsPDDNOS.LTA" = rownames(mout3[which(mout3[,"significance.PDDNOS.NTvsPDDNOS.LTA"] =="DE"),])
)

rapply(listDEGs, length, how="list")

#######################################################################
#for KEGG convert to entreze
#######################################################################

#BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
library(xlsx)

listDEGs2  <- NULL

for (i in 1:length(names(listDEGs))){

  e <- bitr(listDEGs[[i]],fromType='ENSEMBL',toType="ENTREZID", OrgDb="org.Hs.eg.db")
  head(e$ENTREZID)
  
  #save to output list:
  listDEGs2[[i]] <- e$ENTREZID
  
}

names(listDEGs2) <- names(listDEGs)
universe_enterz <-   e <- bitr(GeneSummaryCPM$GeneID,fromType='ENSEMBL',toType="ENTREZID", OrgDb="org.Hs.eg.db")


ck <- compareCluster(geneCluster = listDEGs2,
                     organism="hsa", 
                     universe      = universe_enterz$ENTREZID,
                     fun = "enrichKEGG",
                     qvalueCutoff  = 0.05)




pdf('KeggEnrichment_allDEGs.pdf',w=9,h=11)
d <- dotplot(ck,showCategory = 20)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

ck <- as.data.frame(ck)
ckwrite.csv(ck,file="KEGGenrichments_allDEGs.csv")


#######################################################################
#GO term enrichment
#######################################################################

#######################################################################
#cellular component
#######################################################################

cc <- compareCluster(geneCluster = listDEGs,
                     universe      = GeneSummary2$GeneID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     keyType       = 'ENSEMBL',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     fun = "enrichGO")

#cc <- compareCluster(geneCluster = ck, use_internal_data = TRUE, fun = "enrichKEGG")

  

  


head(cc)
dim(cc)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
cc2 <- simplify(cc, cutoff=0.7, by="p.adjust", select_fun=min)
dim(cc2)


pdf('GO_CCEnrichment_allDEGs.pdf',w=9,h=11)
d <- dotplot(cc2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

cc2 <- as.data.frame(cc2)
write.csv(cc2,file="GO_CCenrichments_allDEGs.csv")

#######################################################################
#biological process
#######################################################################

bp <- compareCluster(geneCluster = listDEGs,
                     universe      = GeneSummary2$GeneID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     keyType       = 'ENSEMBL',
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE,
                     fun = "enrichGO")
head(bp)
dim(bp)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
dim(bp2)


pdf('GO_BPEnrichment_allDEGs.pdf',w=9,h=14)
d <- dotplot(bp2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

bp2 <- as.data.frame(bp2)
write.csv(bp2,file="GO_BPenrichments_allDEGs.csv")

#######################################################################
#molecular function
#######################################################################

mf <- compareCluster(geneCluster = listDEGs,
                     universe      = GeneSummary2$GeneID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "mf",
                     pAdjustMethod = "BH",
                     keyType       = 'ENSEMBL',
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE,
                     fun = "enrichGO")
head(mf)
dim(mf)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
mf2 <- simplify(mf, cutoff=0.7, by="p.adjust", select_fun=min)
dim(mf2)


pdf('GO_MFEnrichment_allDEGs.pdf',w=9,h=14)
d <- dotplot(mf2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

mf2 <- as.data.frame(mf2)
write.csv(mf2,file="GO_MFenrichments_allDEGs.csv")


################################################################################################
#GLIMMA interactive plot building
################################################################################################
#BiocManager::install("Glimma", version = "3.8")
#http://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/Glimma.pdf
library(Glimma)
#setwd("/Users/annieciernia/Box Sync/Ashwood lab collaboration/Experiments/HumanMacrophageLPS/LexogenRNAseq/")

#MDS plots
samplefiles <- v$targets$files

#match ordre of df to v
group.df <- keep.samples[match(samplefiles, keep.samples$filenames),]
group.df <- group.df %>% dplyr::select(Group,Diagnosis,treatment,sex,Age)

levels(group.df$Group)
levels <- c("TD.NT","TD.LPS","TD.LTA",
            "ASD.NT","ASD.LPS", "ASD.LTA",
            "PDDNOS.NT","PDDNOS.LPS", "PDDNOS.LTA")
group.df$Group <- factor(group.df$Group, levels= levels)

#writes out html file:
glMDSPlot(v, groups=group.df)

#DE genes
tfit <- treat(fit2,lfc=log2(2))
results <- decideTests(tfit,adjust.method="BH", method="separate")
summary(results)

#glMDPlot(tfit, counts=v$E,transform=FALSE,status = results,coef=1,
#         anno=tfit$genes, groups = group.df$Group,
#          main=paste("MD plot:",colnames(results)[1]),
 #       display.columns=c("hgnc_symbol", "ensembl_gene_id"),
  #       folder=paste(colnames(results)[1]),side.main="hgnc_symbol")

sample.cols <- c("purple","orange","green","blue","pink",
                 "brown","red","black","grey")[group.df$Group]
for (COEF in 1:6) {
  glMDPlot(tfit, counts=v$E,transform=FALSE,
           anno=tfit$genes, groups = group.df$Group,
           status=results, coef=COEF, main=colnames(tfit)[COEF],
           side.ylab="Log2CPM Expression", side.main="hgnc_symbol",
           folder="glimma_results",sample.cols =sample.cols,
           html = paste("MD-Plot",colnames(results)[COEF]))
}


#comparison 1 vs 2: Typical vs ASD LPS
#genes that are DE in only Typical comparison are highlighted in purple, 

#genes that are DE in both comparisons are highlighted in red.

# #for LPSvsNT Typical vs ASD:
# NTvsASD <- results[,c(3,7)]
# dt2 <- rep(0, nrow(NTvsASD))
# #define rows(genes) where only one conditoin is significant = -1 > purple
# dt2[rowSums(NTvsASD!=0)==1] <- -1
# #define rows(genes) where only both conditions are significant = -1 > red
# dt2[rowSums(NTvsASD!=0)==2] <- 1
# table(dt2)
# colnames(tfit$coef)
# 
# #get genes and write to file
# DFgenes <- cbind(NTvsASD,dt2) 
# Onecondgenes <- DFgenes[which(dt2 == -1),]
# Twocondgenes <- DFgenes[which(dt2 == 1),]
# OneConSig <- rownames(Onecondgenes)
# TwoConSig <- rownames(Twocondgenes)
# 
# DF <- cbind(tfit$coef[,3], y=tfit$coef[,7],tfit$genes)
# colnames(DF)[1:2] <- c(paste("logFC",colnames(tfit$coef)[3], sep=" "),paste("logFC",colnames(tfit$coef)[7], sep=" "))
# 
# DF$Comparison <- "NonDE"
# DF$Comparison[which(DF$ensembl_gene_id %in% OneConSig)] <-c("significant only in one condition")
# DF$Comparison[which(DF$ensembl_gene_id %in% TwoConSig)] <-c("significant in both conditions")
# 
# table(DF$Comparison)
# 
# write.xlsx(DF, "LogFC_LPSvsNT_ASDvsTypical_genes.xlsx")
# 
# cols <- c("purple", "grey", "red")
# glXYPlot(tfit$coef[,3], y=tfit$coef[,7], xlab=paste(colnames(tfit$coef)[3]), 
#          ylab=paste((colnames(tfit$coef)[7])),
#          counts=v$E,transform=FALSE,
#          status=dt2, cols=cols, 
#          anno=tfit$genes, groups = group.df$Group,
#          side.main="hgnc_symbol",side.ylab="Log2CPM Expression", 
#          main="logFCs",
#          folder="glimma_XYplots",sample.cols =sample.cols,
#          html = paste("XY-Plot", colnames(tfit$coef)[3],"vs",colnames(tfit$coef)[7]))
# 
# 
# 
# #for LPSvsNT Typical vs PDDNOS:
# NTvsASD <- results[,c(3,5)]
# dt2 <- rep(0, nrow(NTvsASD))
# #define rows(genes) where only one conditoin is significant = -1 > purple
# dt2[rowSums(NTvsASD!=0)==1] <- -1
# #define rows(genes) where only both conditions are significant = -1 > red
# dt2[rowSums(NTvsASD!=0)==2] <- 1
# table(dt2)
# colnames(tfit$coef)
# 
# #get genes and write to file
# DFgenes <- cbind(NTvsASD,dt2) 
# Onecondgenes <- DFgenes[which(dt2 == -1),]
# Twocondgenes <- DFgenes[which(dt2 == 1),]
# OneConSig <- rownames(Onecondgenes)
# TwoConSig <- rownames(Twocondgenes)
# 
# DF <- cbind(tfit$coef[,3], y=tfit$coef[,5],tfit$genes)
# colnames(DF)[1:2] <- c(paste("logFC",colnames(tfit$coef)[3], sep=" "),paste("logFC",colnames(tfit$coef)[5], sep=" "))
# 
# DF$Comparison <- "NonDE"
# DF$Comparison[which(DF$ensembl_gene_id %in% OneConSig)] <-c("significant only in one condition")
# DF$Comparison[which(DF$ensembl_gene_id %in% TwoConSig)] <-c("significant in both conditions")
# 
# table(DF$Comparison)
# 
# write.xlsx(DF, "LogFC_LPSvsNT_PDDNOSvsTypical_genes.xlsx")
# 
# 
# cols <- c("purple", "grey", "red")
# glXYPlot(tfit$coef[,3], y=tfit$coef[,5], xlab=paste(colnames(tfit$coef)[3]), 
#          ylab=paste((colnames(tfit$coef)[5])),
#          counts=v$E,transform=FALSE,
#          status=dt2, cols=cols, 
#          anno=tfit$genes, groups = group.df$Group,
#          side.main="hgnc_symbol",side.ylab="Log2CPM Expression", 
#          main="logFCs",
#          folder="glimma_XYplots",sample.cols =sample.cols,
#          html = paste("XY-Plot", colnames(tfit$coef)[3],"vs",colnames(tfit$coef)[5]))
# 
# #for LTAvsNT Typical vs ASD:
# NTvsASD <- results[,c(4,8)]
# dt2 <- rep(0, nrow(NTvsASD))
# #define rows(genes) where only one conditoin is significant = -1 > purple
# dt2[rowSums(NTvsASD!=0)==1] <- -1
# #define rows(genes) where only both conditions are significant = -1 > red
# dt2[rowSums(NTvsASD!=0)==2] <- 1
# table(dt2)
# colnames(tfit$coef)
# 
# #get genes and write to file
# DFgenes <- cbind(NTvsASD,dt2) 
# Onecondgenes <- DFgenes[which(dt2 == -1),]
# Twocondgenes <- DFgenes[which(dt2 == 1),]
# OneConSig <- rownames(Onecondgenes)
# TwoConSig <- rownames(Twocondgenes)
# 
# DF <- cbind(tfit$coef[,4], y=tfit$coef[,8],tfit$genes)
# colnames(DF)[1:2] <- c(paste("logFC",colnames(tfit$coef)[4], sep=" "),paste("logFC",colnames(tfit$coef)[8], sep=" "))
# 
# DF$Comparison <- "NonDE"
# DF$Comparison[which(DF$ensembl_gene_id %in% OneConSig)] <-c("significant only in one condition")
# DF$Comparison[which(DF$ensembl_gene_id %in% TwoConSig)] <-c("significant in both conditions")
# 
# table(DF$Comparison)
# 
# write.xlsx(DF, "LogFC_LTAvsNT_ASDvsTypical_genes.xlsx")
# 
# cols <- c("purple", "grey", "red")
# glXYPlot(tfit$coef[,4], y=tfit$coef[,8], xlab=paste(colnames(tfit$coef)[4]), 
#          ylab=paste((colnames(tfit$coef)[8])),
#          counts=v$E,transform=FALSE,
#          status=dt2, cols=cols, 
#          anno=tfit$genes, groups = group.df$Group,
#          side.main="hgnc_symbol",side.ylab="Log2CPM Expression", 
#          main="logFCs",
#          folder="glimma_XYplots",sample.cols =sample.cols,
#          html = paste("XY-Plot", colnames(tfit$coef)[4],"vs",colnames(tfit$coef)[8]))
# 
# #for LTAvsNT Typical vs PDDNOS:
# NTvsASD <- results[,c(4,6)]
# dt2 <- rep(0, nrow(NTvsASD))
# #define rows(genes) where only one conditoin is significant = -1 > purple
# dt2[rowSums(NTvsASD!=0)==1] <- -1
# #define rows(genes) where only both conditions are significant = -1 > red
# dt2[rowSums(NTvsASD!=0)==2] <- 1
# table(dt2)
# colnames(tfit$coef)
# 
# #get genes and write to file
# DFgenes <- cbind(NTvsASD,dt2) 
# Onecondgenes <- DFgenes[which(dt2 == -1),]
# Twocondgenes <- DFgenes[which(dt2 == 1),]
# OneConSig <- rownames(Onecondgenes)
# TwoConSig <- rownames(Twocondgenes)
# 
# DF <- cbind(tfit$coef[,4], y=tfit$coef[,6],tfit$genes)
# colnames(DF)[1:2] <- c(paste("logFC",colnames(tfit$coef)[4], sep=" "),paste("logFC",colnames(tfit$coef)[6], sep=" "))
# 
# DF$Comparison <- "NonDE"
# DF$Comparison[which(DF$ensembl_gene_id %in% OneConSig)] <-c("significant only in one condition")
# DF$Comparison[which(DF$ensembl_gene_id %in% TwoConSig)] <-c("significant in both conditions")
# 
# table(DF$Comparison)
# 
# write.xlsx(DF, "LogFC_LTAvsNT_PDDNOSvsTypical_genes.xlsx")
# 
# cols <- c("purple", "grey", "red")
# glXYPlot(tfit$coef[,4], y=tfit$coef[,6], xlab=paste(colnames(tfit$coef)[4]), 
#          ylab=paste((colnames(tfit$coef)[6])),
#          counts=v$E,transform=FALSE,
#          status=dt2, cols=cols, 
#          anno=tfit$genes, groups = group.df$Group,
#          side.main="hgnc_symbol",side.ylab="Log2CPM Expression", 
#          main="logFCs",
#          folder="glimma_XYplots",sample.cols =sample.cols,
#          html = paste("XY-Plot", colnames(tfit$coef)[4],"vs",colnames(tfit$coef)[6]))
# 


################################################################################################
#Heatmaps
################################################################################################

################################################################################################
#LPS or LTA responsive only in ASD, PDDNOS or Typical
################################################################################################
#remove NA in both LPS and LTA
#mout3 <- read.csv("AllDEG_AllConditions.csv")
#combinedVenn <- mout3[!with(mout3,is.na(mout3$LPSVennGroup)& is.na(mout3$LTAVennGroup)),]

#combinedVenn$combinedVenn <- paste(combinedVenn$LPSVennGroup, combinedVenn$LTAVennGroup, sep = "/")
#remove normally regulated genes in both conditions:
#Typical.NTgreaterLPS, ASD.NTgreaterLPS	Typical.NTgreaterLTA, ASD.NTgreaterLTA
#	Typical.NTlessLPS, ASD.NTlessLPS	Typical.NTlessLTA, ASD.NTlessLTA

#combinedVenn2 <- combinedVenn[!with(combinedVenn,(LPSVennGroup =="Typical.NTgreaterLPS, ASD.NTgreaterLPS")
 #                                  & (LTAVennGroup=="Typical.NTgreaterLTA, ASD.NTgreaterLTA")),]


#combinedVenn2 <- combinedVenn2[!with(combinedVenn2,(LTAVennGroup =="Typical.NTlessLTA, ASD.NTlessLTA")
 #                                   & (LTAVennGroup=="Typical.NTlessLTA, ASD.NTlessLTA")),]

#remove normally regulated genes in LPS condition

#remove NA and subset for heatmap:
combinedVenn <- mout3
combinedVenn$LPSVennGroup <- gsub("less","<",combinedVenn$LPSVennGroup)
combinedVenn$LPSVennGroup <- gsub("greater",">",combinedVenn$LPSVennGroup)
combinedVenn$LTAVennGroup <- gsub("less","<",combinedVenn$LTAVennGroup)
combinedVenn$LTAVennGroup <- gsub("greater",">",combinedVenn$LTAVennGroup)


table(combinedVenn$LTAVennGroup)
write.csv(combinedVenn$LTAVennGroup, "combinedvennLTA.csv")

LPSDEGs <- combinedVenn %>% 
  filter(LPSVennGroup == "PDDNOS.NT>LPS" | LPSVennGroup == "PDDNOS.NT<LPS" | LPSVennGroup == "ASD.NT>LPS" | LPSVennGroup == "ASD.NT<LPS" | LPSVennGroup == "Typical.NT>LPS" | LPSVennGroup == "Typical.NT<LPS") %>%
  arrange(LPSVennGroup)

table(LPSDEGs$LPSVennGroup)
write.csv(combinedVenn$LPSVennGroup, "combinedvennLPS.csv")


LTADEGs <- combinedVenn %>% 
  filter(LTAVennGroup == "PDDNOS.NT>LTA" | LTAVennGroup == "PDDNOS.NT<LTA" | LTAVennGroup == "ASD.NT>LTA" | LTAVennGroup == "ASD.NT<LTA" | LTAVennGroup == "Typical.NT>LTA" | LTAVennGroup == "Typical.NT<LTA") %>%
  arrange(LTAVennGroup)

table(LTADEGs$LTAVennGroup)


################################################################################################
#LPS only genes
################################################################################################


#remove normally regulated DEGs
write.csv(LPSDEGs, file="LPSHeatmapGenes.csv")

#heatmap:
library(pheatmap)

heatDF2 <- LPSDEGs %>%
  dplyr::select(hgnc_symbol, LPSVennGroup, logFC.Typical.NTvsTD.LPS,logFC.ASD.NTvsASD.LPS,logFC.PDDNOS.NTvsPDDNOS.LPS) %>%
  distinct() %>%
  filter(hgnc_symbol != "") %>%
  arrange(LPSVennGroup)



rownames(heatDF2) <- heatDF2$hgnc_symbol

#fix col names

colnames(heatDF2) <- gsub("logFC.","",colnames(heatDF2),fixed=TRUE)
colnames(heatDF2) <- gsub("vsTypical.","vs",colnames(heatDF2),fixed=TRUE)
colnames(heatDF2) <- gsub("vsASD.","vs",colnames(heatDF2),fixed=TRUE)
colnames(heatDF2) <- gsub("vsPDDNOS.","vs",colnames(heatDF2),fixed=TRUE)
names(heatDF2)[3] <- "Typical.NTvsLPS"




#matrix
matrix <- heatDF2 %>% dplyr::select(-LPSVennGroup,-hgnc_symbol)

#define row groups:
annotation_row <- heatDF2[,c(2)]
annotation_row <- as.data.frame(annotation_row)
colnames(annotation_row) <- c("LPS Response")

rownames(annotation_row) = rownames(heatDF2)

head(annotation_row)
#col


annotation_col <- as.data.frame(colnames(matrix))
names(annotation_col) <- c("col")
annotation_col$col <- as.character(annotation_col$col)
annotation_col$col <- gsub("\\."," ",annotation_col$col)

annotation_col <- tidyr::separate(annotation_col,col,into=c("Diagnosis", "Treatment"),sep=" ") 


rownames(annotation_col) = colnames(matrix)

head(annotation_col)

# Specify colors
# #https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
# library(RColorBrewer)
# n <- 27
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# #pie(rep(1,n), col=sample(col_vector, n))

#gradient
#colfunc_blue <- colorRampPalette(c("lightblue", "darkblue"))
#colfunc_green <- colorRampPalette(c("lightgreen", "darkgreen"))

#media,LPS1,LPS2
#cbPalette <- c("royalblue3","darkorange1","brown4") #blue and orange http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# treatment
Var1        <- c("royalblue3","coral1", "orange2")
names(Var1) <- unique(annotation_col$Diagnosis)

#strain
#Var3        <- c("darkseagreen4")
#(Var3) <- unique(annotation_col$Treatment)

Var4 <- c("gold","darkred","blue","darkorange", "darkgreen", "darkorchid4")
names(Var4) <- unique(annotation_row$`LPS Response`)


#Var5 <- c("gold","blue","darkorchid4","grey","darkred","darkorange","darkgreen")
#names(Var5) <- unique(annotation_row$`LTA Response`)


#combined for heatmap labels
anno_colors <- list(Diagnosis=Var1,
                    `LPS Response`=Var4)


#my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

pdf(file = "Heatmap_uniqueLPSDEGs.pdf", wi = 8, he = 12)
pheatmap::pheatmap(matrix, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                    annotation_row = annotation_row,
                   show_rownames = F,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=16, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "Z-score Log2Fold Change Gene Expression")

dev.off()


################################################################################################
#LTA only genes
################################################################################################


#remove normally regulated DEGs
write.csv(LTADEGs, file="LTAHeatmapGenes.csv")

#heatmap:
library(pheatmap)

heatDF2 <- LTADEGs %>%
  dplyr::select(hgnc_symbol, LTAVennGroup, logFC.Typical.NTvsTD.LTA, logFC.ASD.NTvsASD.LTA,logFC.PDDNOS.NTvsPDDNOS.LTA) %>%  
  distinct() %>%
  filter(hgnc_symbol != "") %>%
  arrange(LTAVennGroup)


rownames(heatDF2) <- heatDF2$hgnc_symbol

#fix col names
colnames(heatDF2) <- gsub("logFC.","",colnames(heatDF2),fixed=TRUE)
colnames(heatDF2) <- gsub("vsTypical.","vs",colnames(heatDF2),fixed=TRUE)
colnames(heatDF2) <- gsub("vsASD.","vs",colnames(heatDF2),fixed=TRUE)
colnames(heatDF2) <- gsub("vsPDDNOS.","vs",colnames(heatDF2),fixed=TRUE)
names(heatDF2)[3] <- "Typical.NTvsLTA"



#matrix
matrix <- heatDF2 %>% dplyr::select(-LTAVennGroup,-hgnc_symbol)

#define row groups:
annotation_row <- heatDF2[,c(2)]
annotation_row <- as.data.frame(annotation_row)
colnames(annotation_row) <- c("LTA Response")

rownames(annotation_row) = rownames(heatDF2)

head(annotation_row)
#col


annotation_col <- as.data.frame(colnames(matrix))
names(annotation_col) <- c("col")
annotation_col$col <- as.character(annotation_col$col)
annotation_col$col <- gsub("\\."," ",annotation_col$col)

annotation_col <- tidyr::separate(annotation_col,col,into=c("Diagnosis","Treatment"),sep=" ") 

rownames(annotation_col) = colnames(matrix)

head(annotation_col)

# Specify colors
# #https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
# library(RColorBrewer)
# n <- 27
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# #pie(rep(1,n), col=sample(col_vector, n))

#gradient
#colfunc_blue <- colorRampPalette(c("lightblue", "darkblue"))
#colfunc_green <- colorRampPalette(c("lightgreen", "darkgreen"))

#media,LPS1,LPS2
#cbPalette <- c("royalblue3","darkorange1","brown4") #blue and orange http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

# treatment
Var1        <- c("royalblue3","coral1", "orange2")
names(Var1) <- unique(annotation_col$Diagnosis)

#strain
Var3        <- c("darkseagreen4")
names(Var3) <- unique(annotation_col$Treatment)

Var5 <- c("gold","darkred","blue","darkorange", "darkgreen", "darkorchid4")
names(Var5) <- unique(annotation_row$`LTA Response`)




#combined for heatmap labels
anno_colors <- list(Diagnosis=Var1,
                    Treatment=Var3,
                    `LTA Response`=Var5)


#my_palette <- colorRampPalette(c("blue","white", "red"))(n = 299)

pdf(file = "Heatmap_uniqueLTADEGs.pdf", wi = 8, he = 12)
pheatmap::pheatmap(matrix, 
                   cluster_row = F,
                   cluster_cols = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   show_rownames = F,
                   annotation_names_row=T,
                   #color = my_palette, 
                   fontsize = 16,
                   fontsize_row=16, 
                   fontsize_col = 16,
                   annotation_colors = anno_colors, 
                   scale = c("row"),
                   main = "LogFold Change Gene Expression")

dev.off()



################################################################################################
################################################################################################
################################################################################################
#Enrichments LPS directional
################################################################################################
################################################################################################
library(clusterProfiler)
#The input for geneCluster parameter should be a named list of gene IDs. 
LPSDEGs$GeneID <- as.character(LPSDEGs$GeneID)
LTADEGs$GeneID <- as.character(LTADEGs$GeneID)


#list of ensembl ids for GO
LPS_ensembl <- LPSDEGs %>% dplyr::select(LPSVennGroup,GeneID) %>% na.omit
LTA_ensembl <- LTADEGs %>% dplyr::select(LTAVennGroup,GeneID) %>% na.omit

#list of entrezid for KEGG
LPS_entrez <- subset(mout3, GeneID %in% LPSDEGs$GeneID )
LPS_entrez <- LPS_entrez %>% dplyr::select(LPSVennGroup,entrezgene_id) %>% na.omit

LTA_entrez <- subset(mout3, GeneID %in% LTADEGs$GeneID )
LTA_entrez <- LTA_entrez %>% dplyr::select(LTAVennGroup,entrezgene_id) %>% na.omit

#fix names
LPS_entrez$LPSVennGroup <- gsub("greater"," > ",LPS_entrez$LPSVennGroup)
LPS_entrez$LPSVennGroup <- gsub("less"," < ",LPS_entrez$LPSVennGroup)

LTA_entrez$LTAVennGroup <- gsub("greater"," > ",LTA_entrez$LTAVennGroup)
LTA_entrez$LTAVennGroup <- gsub("less"," < ",LTA_entrez$LTAVennGroup)


  
#convert each data frome to a list of lists:
 #  Split on LPSVennGroup
LPS_entrez_list <- split(LPS_entrez$entrezgene_id, f = LPS_entrez$LPSVennGroup )
LTA_entrez_list <- split(LTA_entrez$entrezgene_id, f = LTA_entrez$LTAVennGroup )
LPS_ensembl_list <- split(LPS_ensembl$GeneID, f = LPS_ensembl$LPSVennGroup )
LTA_ensembl_list <- split(LTA_ensembl$GeneID, f = LTA_ensembl$LTAVennGroup )


#master lists
#list_entrez <- append(LPS_entrez_list,LTA_entrez_list)
#list_ensembl <- append(LPS_ensembl_list,LTA_ensembl_list)

#number of genes perlist:
count <- rapply(LPS_ensembl_list, length, how="list")
write.csv(count,"NumberGenes_LPSHeatmapGroups.csv")

count <- rapply(LTA_ensembl_list, length, how="list")
write.csv(count,"NumberGenes_LTAHeatmapGroups.csv")


#######################################################################
#for KEGG convert to entreze
#######################################################################
#Background:
universe_enterz <- mout3 %>% dplyr::select(entrezgene_id) %>% na.omit
universe_enterz$entrezgene_id <- as.character(universe_enterz$entrezgene_id) #11418

#LPS
kegg_LPS <- compareCluster(geneCluster = LPS_entrez_list,
                           organism="hsa", 
                           universe      = universe_enterz$entrezgene_id,
                           fun = "enrichKEGG",
                           qvalueCutoff  = 0.05
                         )


pdf('LPSHeatmap_KeggEnrichment.pdf',w=8,h=4)
d <- dotplot(kegg_LPS,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

kegg_LPS <- as.data.frame(kegg_LPS)
write.csv(kegg_LPS,file="LPSHeatmapKEGGenrichments.csv")


#LTA
kegg_LTA <- compareCluster(geneCluster = LTA_entrez_list,
                           organism="hsa", 
                           universe      = universe_enterz$entrezgene_id,
                           fun = "enrichKEGG",
                           qvalueCutoff  = 0.05)

pdf('LTAHeatmap_KeggEnrichment.pdf',w=8,h=4)
d <- dotplot(kegg_LTA,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

kegg_LTA <- as.data.frame(kegg_LTA)
write.csv(kegg_LTA,file="LTAHeatmapKEGGenrichments.csv")


#######################################################################
#GO term enrichment
#######################################################################
#Background:
universe_ensembl <- mout3 %>% dplyr::select(GeneID) %>% na.omit
universe_ensembl$GeneID <- as.character(universe_ensembl$GeneID)

#######################################################################
#cellular component
#######################################################################
library(clusterProfiler)
#LPS
LPScc <- compareCluster(geneCluster = LPS_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LPScc)
dim(LPScc)

#didn't simplify due to low enrichment
#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
#LPScc2 <- clusterProfiler::simplify(LPScc, cutoff=0.7, by="p.adjust", select_fun=min)
#dim(LPScc2)


pdf('LPSHeatmap_GO_CCEnrichment.pdf',w=8,h=5)
d <- dotplot(LPScc,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LPScc <- as.data.frame(LPScc)
write.csv(LPScc2,file="LPSHeatmap_GO_CCenrichments.csv")

#LTA
LTAcc <- compareCluster(geneCluster = LTA_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LTAcc)
dim(LTAcc)

#didn't simplify due to low enrichment
#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
#LTAcc2 <- clusterProfiler::simplify(LTAcc, cutoff=0.7, by="p.adjust", select_fun=min)
#dim(LTAcc2)


pdf('LTAHeatmap_GO_CCEnrichment.pdf',w=8,h=5)
d <- dotplot(LTAcc,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LTAcc <- as.data.frame(LTAcc)
write.csv(LTAcc,file="LTAHeatmap_GO_CCenrichments.csv")

#######################################################################
#biological process
#######################################################################
LPSBP <- compareCluster(geneCluster = LPS_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.1,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LPSBP)
dim(LPSBP)

#didn't simplify due to low enrichment
#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
#LPSBP2 <- clusterProfiler::simplify(LPSBP, cutoff=0.7, by="p.adjust", select_fun=min)
#dim(LPSBP2)


pdf('LPSHeatmap_GO_BPEnrichment.pdf',w=8,h=5)
d <- dotplot(LPSBP2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LPSBP <- as.data.frame(LPSBP)
write.csv(LPSBP,file="LPSHeatmap_GO_BPenrichments.csv")

#LTA
LTABP <- compareCluster(geneCluster = LTA_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LTABP)
dim(LTABP)

#didn't simplify due to low enrichment
#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
#LTABP2 <- clusterProfiler::simplify(LTABP, cutoff=0.7, by="p.adjust", select_fun=min)
#dim(LTABP2)


pdf('LTAHeatmap_GO_BPEnrichment.pdf',w=8,h=5)
d <- dotplot(LTABP,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LTABP2 <- as.data.frame(LTABP)
write.csv(LTABP,file="LTAHeatmap_GO_BPenrichments.csv")

#######################################################################
#molecular function
#######################################################################
LPSMF <- compareCluster(geneCluster = LPS_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = enrichGO)
head(LPSMF)
dim(LPSMF)

#No enrichment found in any of gene cluster
#didn't simplify due to low enrichment
#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
#LPSMF2 <- clusterProfiler::simplify(LPSMF, cutoff=0.7, by="p.adjust", select_fun=min)
#dim(LPSMF2)


#pdf('LPSHeatmap_GO_MFEnrichment.pdf',w=8,h=5)
#d <- dotplot(LPSMF,showCategory = 40)
#d + theme(axis.text.x = element_text(angle = 45, hjust=1))
#dev.off()

#LPSMF2 <- as.data.frame(LPSMF)
#write.csv(LPSMF,file="LPSHeatmap_GO_MFenrichments.csv")

#LTA
LTAMF <- compareCluster(geneCluster = LTA_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LTAMF)
dim(LTAMF)

#No enrichment found in any of gene cluster
#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
# LTAMF2 <- clusterProfiler::simplify(LTAMF, cutoff=0.7, by="p.adjust", select_fun=min)
# dim(LPSMF2)

# pdf('LTAHeatmap_GO_MFEnrichment.pdf',w=8,h=5)
# d <- dotplot(LTAMF2,showCategory = 40)
# d + theme(axis.text.x = element_text(angle = 45, hjust=1))
# dev.off()
# 
# LTAMF2 <- as.data.frame(LPSMF2)
# write.csv(LTAMF2,file="LTAHeatmap_GO_MFenrichments.csv")
##################################################################
#Upset plots uniquely regulated genes
##################################################################

library("UpSetR")

#list of conditions to overlap
listInput = append(LPS_ensembl_list,LTA_ensembl_list)

# listInput = append(x,x2)
# rapply(listInput, length, how="list")
# 
# names(listInput) <- gsub("greater",">",names(listInput))
# names(listInput) <- gsub("less","<",names(listInput))

#colors:
queries = 
  list( #only in NT
    list(query =  intersects, params = list("Typical.NT<LPS"), color = "blue", active = T,query.name = c("NT Increase only")), 
    list(query =  intersects, params = list("Typical.NT<LTA"), color = "blue", active = T,query.name = c("NT Increase only")), 
    list(query =  intersects, params = list("Typical.NT<LPS", "Typical.NT<LTA"), color = "blue", active = T,query.name = c("NT Increase only")), 
    
    list(query =  intersects, params = list("Typical.NT>LPS"), color = "darkorange", active = T,query.name = c("NT Decrease only")), 
    list(query =  intersects, params = list("Typical.NT>LTA"), color = "darkorange", active = T,query.name = c("NT Decrease only")), 
    list(query =  intersects, params = list("Typical.NT>LPS", "Typical.NT>LTA"), color = "darkorange", active = T,query.name = c("NT Decrease only")), 
    
    
    
    #only in ASD LPS or LTA
    list(query =  intersects, params = list("ASD.NT<LPS"), color = "purple", active = T,query.name = c("NT Increase only")), 
    list(query =  intersects, params = list("ASD.NT<LTA"), color = "purple", active = T,query.name = c("NT Increase only")), 
    list(query =  intersects, params = list("ASD.NT<LPS", "ASD.NT<LTA"), color = "purple", active = T,query.name = c("NT Increase only")), 
    
    list(query =  intersects, params = list("ASD.NT>LPS"), color = "yellow", active = T,query.name = c("NT Decrease only")), 
    list(query =  intersects, params = list("ASD.NT>LTA"), color = "yellow", active = T,query.name = c("NT Decrease only")), 
    list(query =  intersects, params = list("ASD.NT>LPS", "ASD.NT>LTA"), color = "yellow", active = T,query.name = c("NT Decrease only")),
    
    #only in PDDNOS LPS or LTA
    list(query =  intersects, params = list("PDDNOS.NT<LPS"), color = "gold", active = T,query.name = c("NT Increase only")), 
    list(query =  intersects, params = list("PDDNOS.NT<LTA"), color = "gold", active = T,query.name = c("NT Increase only")), 
    list(query =  intersects, params = list("PDDNOS.NT<LPS", "PDDNOS.NT<LTA"), color = "gold", active = T,query.name = c("NT Increase only")), 
    
    list(query =  intersects, params = list("PDDNOS.NT>LPS"), color = "darkred", active = T,query.name = c("NT Decrease only")), 
    list(query =  intersects, params = list("PDDNOS.NT>LTA"), color = "darkred", active = T,query.name = c("NT Decrease only")), 
    list(query =  intersects, params = list("PDDNOS.NT>LPS", "PDDNOS.NT>LTA"), color = "darkred", active = T,query.name = c("NT Decrease only"))
    
  )

#blue = genes regulated only in PDDNOS
#red = genes regulated only in ASD
#orange = genes regulated only in NT
#magenta = genes only in ASD or PDDNOS

pdf(file = "UpSetPlot_uniqueDEGsgreaterthan.pdf", wi = 8, he = 7, useDingbats=F)
upset(fromList(listInput),sets = c("PDDNOS.NT>LPS", "PDDNOS.NT>LTA", "ASD.NT>LTA", "ASD.NT>LPS", "Typical.NT>LPS", "Typical.NT>LTA"), 
      order.by = "freq",
      #number.angles =10,
      keep.order = TRUE,
      text.scale = c(2)
     )

dev.off()

pdf(file = "UpSetPlot_uniqueDEGslessthan.pdf", wi = 8, he = 7, useDingbats=F)
upset(fromList(listInput),sets = c("PDDNOS.NT<LPS", "PDDNOS.NT<LTA", "ASD.NT<LTA", "ASD.NT<LPS", "Typical.NT<LPS", "Typical.NT<LTA"), 
      order.by = "freq",
      #number.angles =10,
      keep.order = TRUE,
      text.scale = c(2)
)

dev.off()

#pdf(file = "UpSetPlot_uniqueDEGstypical.pdf", wi = 8, he = 7, useDingbats=F)
#upset(fromList(listInput),sets = c("Typical.NT>LTA", "Typical.NT>LPS", "Typical.NT<LPS", "Typical.NT<LTA"), 
      #order.by = "freq",
      #number.angles =10,
      #keep.order = TRUE,
      #text.scale = c(2)
#)

#dev.off()

#pdf(file = "UpSetPlot_uniqueDEGsASD.pdf", wi = 8, he = 7, useDingbats=F)
#upset(fromList(listInput),sets = c("ASD.NT<LTA", "ASD.NT<LPS", "ASD.NT>LPS", "ASD.NT>LTA"), 
      #order.by = "freq",
      #number.angles =10,
      #keep.order = TRUE,
      #text.scale = c(2)
#)

#dev.off()

#pdf(file = "UpSetPlot_uniqueDEGsPDDNOS.pdf", wi = 8, he = 7, useDingbats=F)
#upset(fromList(listInput),sets = c("PDDNOS.NT<LTA", "PDDNOS.NT<LPS", "PDDNOS.NT>LPS"), 
      #order.by = "freq",
      #number.angles =10,
      #keep.order = TRUE,
      #text.scale = c(2)
#)

#dev.off()

#pdf(file = "UpSetPlot_uniqueDEGs.pdf", wi = 8, he = 7, useDingbats=F)
#upset(fromList(listInput),sets = rev(names(listInput)), 
      #order.by = "freq",
      #number.angles =10,
      #keep.order = TRUE,
     #text.scale = c(2)
#)

#dev.off()

#table(mout3$LPSVennGroup,mout3$LTAVennGroup)

################################################################################################
#Typical NT< LPS vs Typical NT<LTA
################################################################################################

#######################################################################
#DE plots for Upregulated genes with LTA treatment
#######################################################################
# DecrlistTyp <- list(`Typical.NT>LPS` = listInput$`Typical.NT>LPS`,
#                  `Typical.NT>LTA`= listInput$`Typical.NT>LTA`)
# 
# pdf("Venn_UniqueTypical_NTgreaterLPS.NTgreaterLTA.pdf",height=10,width = 10)
# venn(DecrlistTyp)
# dev.off()
# 
# IncrlistTyp <- list(`Typical.NT<LPS` = listInput$`Typical.NT<LPS`,
#                  `Typical.NT<LTA`= listInput$`Typical.NT<LTA`)
# 
# pdf("Venn_UniqueTypical_NTlessLPS.NTlessLTA.pdf",height=10,width = 10)
# venn(IncrlistTyp)
# dev.off()
# 
# DecrlistASD<- list(`ASD.NT>LPS` = listInput$`ASD.NT>LPS`,
#                    `ASD.NT>LTA`= listInput$`ASD.NT>LTA`)
# 
# pdf("Venn_UniqueASD_NTgreaterLPS.NTgreaterLTA.pdf",height=10,width = 10)
# venn(DecrlistASD)
# dev.off()
# 
# IncrlistASD<- list(`ASD.NT<LPS` = listInput$`ASD.NT<LPS`,
#                    `ASD.NT<LTA`= listInput$`ASD.NT<LTA`)
# 
# pdf("Venn_UniqueASD_NTlessLPS.NTlessLTA.pdf",height=10,width = 10)
# venn(IncrlistASD)
# dev.off()
# 
# DecrlistPDDNOS<- list(`PDDNOS.NT>LPS` = listInput$`PDDNOS.NT>LPS`,
#                    `PDDNOS.NT>LTA`= listInput$`PDDNOS.NT>LTA`)
# 
# pdf("Venn_UniquePDDNOS_NTgreaterLPS.NTgreaterLTA.pdf",height=10,width = 10)
# venn(DecrlistPDDNOS)
# dev.off()
# 
# IncrlistPDDNOS<- list(`PDDNOS.NT<LPS` = listInput$`PDDNOS.NT<LPS`,
#                    `PDDNOS.NT<LTA`= listInput$`PDDNOS.NT<LTA`)
# 
# pdf("Venn_UniquePDDNOS_NTlessLPS.NTlessLTA.pdf",height=10,width = 10)
# venn(IncrlistPDDNOS)
# dev.off()
# 
# rapply(DecrlistTyp, length, how="list")
##################################################################
#Fisher's exact test for list overlaps LTA
##################################################################
library(GeneOverlap)
#https://bioc.ism.ac.jp/packages/3.4/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf

#background = all genes expressed in the experiment
BG <- length(GeneSummaryCPM$GeneID)

genelists <- list(IncrlistTyp,DecrlistTyp,IncrlistASD,DecrlistASD,IncrlistPDDNOS, DecrlistPDDNOS)
names(genelists) <- c("IncrlistTyp","DecrlistTyp","IncrlistASD","DecrlistASD","IncrlistPDDNOS","DecrlistPDDNOS")

output <- NULL
for (i in 1:length(genelists)){ 
  print(i)
#overlap list all vs all
go.obj <-  newGOM(genelists[[i]],genome.size=BG)
print(go.obj)

#get p values
options(digits=10)
pval <- getMatrix(go.obj, name="pval")

pval <- as.data.frame(pval)

#get odds ratio
OR <- getMatrix(go.obj, "odds.ratio")
OR <- as.data.frame(OR)

Overlaps <- cbind(pval,OR)
Overlaps$comparison <- rownames(Overlaps)

output <- rbind(Overlaps,output)

}

output$FDR <- p.adjust(output$pval,method="fdr")

write.csv(output,"UniqueDEGsOverlapsFishersExactTest.csv")

#intersecting genes:
#typical down
intersection <- intersect(DecrlistTyp$`Typical.NT>LPS`,DecrlistTyp$`Typical.NT>LTA`)

inter.genes <- mout3[which(mout3$GeneID %in% intersection),]

write.csv(inter.genes,"Genes_TypicalDownOverlapsFishersExactTest.csv")

#typical up
intersection <- intersect(IncrlistTyp$`Typical.NT<LPS`,IncrlistTyp$`Typical.NT<LTA`)

inter.genes <- mout3[which(mout3$GeneID %in% intersection),]
dim(inter.genes)

write.csv(inter.genes,"Genes_TypicalUPOverlapsFishersExactTest.csv")

#intersecting genes:
#ASD down
intersection <- intersect(DecrlistASD$`ASD.NT>LPS`,DecrlistASD$`ASD.NT>LTA`)

inter.genes <- mout3[which(mout3$GeneID %in% intersection),]
dim(inter.genes)
write.csv(inter.genes,"Genes_ASDDownOverlapsFishersExactTest.csv")

#ASD up
intersection <- intersect(IncrlistASD$`ASD.NT<LPS`,IncrlistASD$`ASD.NT<LTA`)

inter.genes <- mout3[which(mout3$GeneID %in% intersection),]
dim(inter.genes)

write.csv(inter.genes,"Genes_ASDUPOverlapsFishersExactTest.csv")

#intersecting genes:
#PDDNOS down
intersection <- intersect(DecrlistASD$`PDDNOS.NT>LPS`,DecrlistPDDNOS$`PDDNOS.NT>LTA`)

inter.genes <- mout3[which(mout3$GeneID %in% intersection),]
dim(inter.genes)
write.csv(inter.genes,"Genes_PDDNOSDownOverlapsFishersExactTest.csv")

#PDDNOS up
intersection <- intersect(IncrlistPDDNOS$`PDDNOS.NT<LPS`,IncrlistPDDNOS$`PDDNOS.NT<LTA`)

inter.genes <- mout3[which(mout3$GeneID %in% intersection),]
dim(inter.genes)

write.csv(inter.genes,"Genes_PDDNOSUPOverlapsFishersExactTest.csv")

################################################################################################
################################################################################################
#repeat enrichments collapsed by direction
################################################################################################
################################################################################################
library(clusterProfiler)
#The input for geneCluster parameter should be a named list of gene IDs. 
LPSDEGs$GeneID <- as.character(LPSDEGs$GeneID)
LTADEGs$GeneID <- as.character(LTADEGs$GeneID)


#list of ensembl ids for GO
LPS_ensembl <- LPSDEGs %>% dplyr::select(LPSVennGroup,GeneID) %>% na.omit
LPS_ensembl$GeneID <- as.character(LPS_ensembl$GeneID)
LTA_ensembl <- LTADEGs %>% dplyr::select(LTAVennGroup,GeneID) %>% na.omit
LTA_ensembl$GeneID <- as.character(LTA_ensembl$GeneID)

#list of entrezid for KEGG
LPS_entrez <- subset(mout3, GeneID %in% LPSDEGs$GeneID )
LPS_entrez <- LPS_entrez %>% dplyr::select(LPSVennGroup,entrezgene_id) %>% na.omit
LPS_entrez$entrezgene_id <- as.character(LPS_entrez$entrezgene_id)


LTA_entrez <- subset(mout3, GeneID %in% LTADEGs$GeneID )
LTA_entrez <- LTA_entrez %>% dplyr::select(LTAVennGroup,entrezgene_id) %>% na.omit
LTA_entrez$entrezgene_id <- as.character(LTA_entrez$entrezgene_id)


#fix names
LPS_entrez$LPSVennGroup <- gsub("greater","vs",LPS_entrez$LPSVennGroup)
LPS_entrez$LPSVennGroup <- gsub("less","vs",LPS_entrez$LPSVennGroup)

LTA_entrez$LTAVennGroup <- gsub("greater","vs",LTA_entrez$LTAVennGroup)
LTA_entrez$LTAVennGroup <- gsub("less","vs",LTA_entrez$LTAVennGroup)

#fix names
LPS_ensembl$LPSVennGroup <- gsub(">","vs",LPS_ensembl$LPSVennGroup)
LPS_ensembl$LPSVennGroup <- gsub("<","vs",LPS_ensembl$LPSVennGroup)

LTA_ensembl$LTAVennGroup <- gsub(">","vs",LTA_ensembl$LTAVennGroup)
LTA_ensembl$LTAVennGroup <- gsub("<","vs",LTA_ensembl$LTAVennGroup)
  
#convert each data frome to a list of lists:
#Split on LPSVennGroup
LPS_entrez_list <- split(LPS_entrez$entrezgene_id, f = LPS_entrez$LPSVennGroup )
LTA_entrez_list <- split(LTA_entrez$entrezgene_id, f = LTA_entrez$LTAVennGroup )
LPS_ensembl_list <- split(LPS_ensembl$GeneID, f = LPS_ensembl$LPSVennGroup )
LTA_ensembl_list <- split(LTA_ensembl$GeneID, f = LTA_ensembl$LTAVennGroup )


#master lists
list_entrez <- append(LPS_entrez_list,LTA_entrez_list)
list_ensembl <- append(LPS_ensembl_list,LTA_ensembl_list)

#number of genes perlist:
count <- rapply(list_entrez, length, how="list")
#write.xlsx(count,"NumberGenes_LPSHeatmapGroups.xlsx")

count <- rapply(list_ensembl, length, how="list")
#write.xlsx(count,"NumberGenes_LTAHeatmapGroups.xlsx")


#######################################################################
#for KEGG convert to entreze
#######################################################################
#Background:
universe_enterz <- mout3 %>% dplyr::select(entrezgene_id) %>% na.omit
universe_enterz$entrezgene_id <- as.character(universe_enterz$entrezgene_id)

#LPS
#LPS
kegg_LPS <- compareCluster(geneCluster = LPS_entrez_list,
                           organism="hsa", 
                           universe      = universe_enterz$entrezgene_id,
                           fun = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05
)


pdf('LPSHeatmap_KeggEnrichment.pdf',w=8,h=4)
d <- dotplot(kegg_LPS,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

kegg_LPS <- as.data.frame(kegg_LPS)
write.csv(kegg_LPS,file="LPSHeatmapKEGGenrichments.csv")


#LTA
kegg_LTA <- compareCluster(geneCluster = LTA_entrez_list,
                           organism="hsa", 
                           universe      = universe_enterz$entrezgene_id,
                           fun = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05)

pdf('LTAHeatmap_KeggEnrichment.pdf',w=8,h=4)
d <- dotplot(kegg_LTA,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

kegg_LTA <- as.data.frame(kegg_LTA)
write.csv(kegg_LTA,file="LTAHeatmapKEGGenrichments.csv")


#######################################################################
#GO term enrichment - changed p-value cut off from 0.01 to 0.05 (7-15-21)
#######################################################################
#Background:
universe_ensembl <- mout3 %>% dplyr::select(GeneID) %>% na.omit
universe_ensembl$GeneID <- as.character(universe_ensembl$GeneID)

#######################################################################
#cellular component
#######################################################################
library(clusterProfiler)
#LPS
LPScc <- compareCluster(geneCluster = LPS_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LPScc)
dim(LPScc)


#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
LPScc2 <- clusterProfiler::simplify(LPScc, cutoff=0.7, by="p.adjust", select_fun=min)
dim(LPScc2)


pdf('LPSHeatmap_GO_CCEnrichment.pdf',w=8,h=5)
d <- dotplot(LPScc2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LPScc2 <- as.data.frame(LPScc2)
write.csv(LPScc2,file="LPSHeatmap_GO_CCenrichments.csv")

#LTA
LTAcc <- compareCluster(geneCluster = LTA_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "CC",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LTAcc)
dim(LTAcc)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
LTAcc2 <- clusterProfiler::simplify(LTAcc, cutoff=0.7, by="p.adjust", select_fun=min)
dim(LTAcc2)


pdf('LTAHeatmap_GO_CCEnrichment.pdf',w=8,h=5)
d <- dotplot(LTAcc2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LTAcc2 <- as.data.frame(LTAcc2)
write.csv(LTAcc2,file="LTAHeatmap_GO_CCenrichments.csv")

#######################################################################
#biological process
#######################################################################
LPSBP <- compareCluster(geneCluster = LPS_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LPSBP)
dim(LPSBP)


#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
LPSBP2 <- clusterProfiler::simplify(LPSBP, cutoff=0.7, by="p.adjust", select_fun=min)
dim(LPSBP2)


pdf('LPSHeatmap_GO_BPEnrichment.pdf',w=8,h=5)
d <- dotplot(LPSBP2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LPSBP2 <- as.data.frame(LPSBP2)
write.csv(LPSBP2,file="LPSHeatmap_GO_BPenrichments.csv")

#LTA
LTABP <- compareCluster(geneCluster = LTA_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LTABP)
dim(LTABP)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
LTABP2 <- clusterProfiler::simplify(LTABP, cutoff=0.7, by="p.adjust", select_fun=min)
dim(LTABP2)


pdf('LTAHeatmap_GO_BPEnrichment.pdf',w=8,h=5)
d <- dotplot(LTABP2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LTABP2 <- as.data.frame(LTABP2)
write.csv(LTABP2,file="LTAHeatmap_GO_BPenrichments.csv")

#######################################################################
#molecular function
#######################################################################
LPSMF <- compareCluster(geneCluster = LPS_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = enrichGO)
head(LPSMF)
dim(LPSMF)

#No enrichment found in any of gene cluster
#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
#LPSMF2 <- clusterProfiler::simplify(LPSMF, cutoff=0.7, by="p.adjust", select_fun=min)
#dim(LPSMF2)


#pdf('LPSHeatmap_GO_MFEnrichment.pdf',w=8,h=5)
#d <- dotplot(LPSMF,showCategory = 40)
#d + theme(axis.text.x = element_text(angle = 45, hjust=1))
#dev.off()

#LPSMF2 <- as.data.frame(LPSMF)
#write.csv(LPSMF,file="LPSHeatmap_GO_MFenrichments.csv")

#LTA
LTAMF <- compareCluster(geneCluster = LTA_ensembl_list,
                        universe      = universe_ensembl$GeneID,
                        OrgDb         = org.Hs.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        keyType       = 'ENSEMBL',
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        fun = "enrichGO")
head(LTAMF)
dim(LTAMF)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
LTAMF2 <- clusterProfiler::simplify(LTAMF, cutoff=0.7, by="p.adjust", select_fun=min)
dim(LTAMF2)

pdf('LTAHeatmap_GO_MFEnrichment.pdf',w=8,h=5)
d <- dotplot(LTAMF2,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

LTAMF2 <- as.data.frame(LTAMF2)
write.csv(LTAMF2,file="LTAHeatmap_GO_MFenrichments.csv")

#Graphs from overlap analysis
#author: Annie Vogel Ciernia
#a.ciernia@gmail.com
#updated 6-21-21 by Megan Rowland
##############################################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(lsmeans)
library(nlme)
library(cowplot)
theme_set(theme_cowplot())
library(xlsx)
##############################################################################################################
#get human gene lists
##############################################################################################################

#get asd gene lists
setwd("C:/Users/Megan/Documents/Ciernia Lab/Bioninformatics/HumanMacrophageAshwood/6-11-21/genelistoverlaps")
ASDoverlaplist <- read.csv("TableS8_MasterOverlapsDMRgenesV3.csv")
unique(ASDoverlaplist$List)
#####################################################
#list info summary
infolist <- ASDoverlaplist %>% dplyr::select(List, List.description,Study,tissue, ensembl_gene_id) %>%
  distinct() 

#remove some lists:
infolist <- infolist[!grepl("isoModules", infolist$List),]
infolist <- infolist[!grepl("Li2018", infolist$List),]
infolist <- infolist[!grepl("geneModules", infolist$List),]
infolist <- infolist[!grepl("NCrnas", infolist$List),]
infolist <- infolist[!grepl("DTU", infolist$List),]
infolist <- infolist[!grepl("DTE", infolist$List),]
infolist <- infolist[!grepl("ASD cortex differentially spliced genes", infolist$List),]
infolist <- infolist[!grepl("Nardone ASD DMR genes", infolist$List),]



unique(infolist$List)

#fix names:
infolist$List <- gsub(" DGE","",infolist$List)
infolist$List <- gsub("genes with LGD mutations in ASD","ASD LGD mutations",infolist$List)
infolist$List <- gsub("LGD Recurrent","ASD LGD recurrent mutations",infolist$List)
infolist$List <- gsub("SFARI","SFARI ASD",infolist$List)
infolist$List <- gsub("Gandal ","", infolist$List)
infolist$List <- gsub("Miller ","", infolist$List)


#combine ID lists:
infolist$List <- gsub("Candidate ID","ID genetic risk",infolist$List)
infolist$List <- gsub("de novo CNV ID","ID genetic risk",infolist$List)
infolist$List <- gsub("ID genes","ID genetic risk",infolist$List)
infolist$List <- gsub("Known ID","ID genetic risk",infolist$List)
infolist$List <- gsub("genes with LGD mutations in ID","ID genetic risk",infolist$List)

unique(infolist$List)

#final name fix
infolist$List <- gsub(" genes","", infolist$List)

#make lists for analysis:

library("GeneOverlap")
ASDoverlaplist2 <- infolist %>% dplyr::select(List,ensembl_gene_id) %>% distinct()
ASDoverlaplist2$ensembl_gene_id <- as.character(ASDoverlaplist2$ensembl_gene_id)
ASDoverlaplist2 <- as.data.frame(ASDoverlaplist2)
List <- split(ASDoverlaplist2$ensembl_gene_id, ASDoverlaplist2$List)
names(List)

#write out summary file for references
infolist2 <- infolist %>% dplyr::select(List, List.description,Study,tissue) %>%
  distinct() 
write.csv(infolist2,file="NDDgenelists_Info.csv")

#background genes: all human genes in hg38
human_genes <- read.csv("Allhg38_genes.csv")

hg38genome <- length(unique(human_genes$ensembl_gene_id))
#64914

##########################################################################################################
# PBMCS DEGS from heatmap
##########################################################################################################
LPSDEGs <- read.csv("C:/Users/Megan/Documents/Ciernia Lab/Bioninformatics/HumanMacrophageAshwood/test/LPSHeatmapGenes.csv")

#LPS groups
LPSDEGs$LPSVennGroup <- as.character(LPSDEGs$LPSVennGroup)
unique(LPSDEGs$LPSVennGroup)

#list info summary
LPSDEGinfo <- LPSDEGs %>% dplyr::select(LPSVennGroup, GeneID) %>%
  distinct() 

colnames(LPSDEGinfo) <- c("VennGroup","ensemble_gene_id")


#LTA groups
LTADEGs <- read.csv("C:/Users/Megan/Documents/Ciernia Lab/Bioninformatics/HumanMacrophageAshwood/test/LTAHeatmapGenes.csv")

#LPS groups
LTADEGs$LTAVennGroup <- as.character(LTADEGs$LTAVennGroup)
unique(LTADEGs$LTAVennGroup)

#list info summary
LTADEGinfo <- LTADEGs %>% dplyr::select(LTAVennGroup, GeneID) %>%
  distinct() 

colnames(LTADEGinfo) <- c("VennGroup","ensemble_gene_id")


DEGs <- rbind(LPSDEGinfo, LTADEGinfo)
DEGs$ensemble_gene_id <- as.character(DEGs$ensemble_gene_id)

#combine to list:
asdlist <- split(DEGs$ensemble_gene_id, DEGs$VennGroup)


##########################################################################################################
#combine lists
##########################################################################################################
Masterlist <- append(List,asdlist)
names(Masterlist)

##########################################################################################################
#select wanted gene lists
##########################################################################################################
library(rlist)

#Masterlist2 <-list.remove(Masterlist, names(Masterlist)[grepl("isoModules", names(Masterlist))])
#names(Masterlist2)



#number of genes in each list:
genecount <- rapply(Masterlist, length, how="list")
genecount <- as.data.frame(unlist(genecount))
genecount



##########################################################################################################
#Overlap Function: Fisher's one tailed exact test
##########################################################################################################
library("GeneOverlap")

#for multiple lists:
genelists <- Masterlist
genomesize <- hg38genome


Overlap_fxn <- function(targetlistname,genelists,genomesize){
  out <- NULL
  target <- genelists[[targetlistname]]
  inputlistname <- names(genelists[targetlistname])
  
  for (i in 1:length(genelists)) { 
    
    #call gene overlaps
    go.obj <- newGeneOverlap(target,
                             genelists[[i]],
                             genome.size=genomesize)
    
    #perform test
    go.obj <- testGeneOverlap(go.obj) #returns onetailed pvalue
    #return odds ratio:
    OR <- getOddsRatio(go.obj)
    pvalue <- getPval(go.obj)
    
    #extract contingency table
    CT <- getContbl(go.obj)
    notAnotB <- CT[1,1]
    inAnotB <- CT[1,2]
    inBnotA <- CT[2,1]
    inBinA <- CT[2,2]
    
    CTlist <- cbind(notAnotB,inAnotB,inBnotA,inBinA)
    
    
    #two sided fisher's exact test
    #test <- fisher.test(CT,alternative='two.sided')
    
    #get gene list B:
    intersection <- go.obj@intersection
    intersection_ensembl <- paste(as.character(intersection),collapse=", ",sep="")
    
    #get intersection gene names
    intersection_genenames <- human_genes$hgnc_symbol[which(human_genes$ensembl_gene_id %in% intersection)]
    intersection_genenames <- paste(as.character(intersection_genenames),collapse=", ",sep="")
    
    #get listname
    listname <- paste(names(genelists[i]))
    
    results <- cbind(listname,pvalue,OR, CTlist,intersection_ensembl,intersection_genenames )
    
    names(results) <- c("listname","pvalue","OR","notAnotB","inAnotB","inBnotA","inBinA","ensembl","geneID")
    out <- rbind(out,results) 
    
  }
  
  #remove first row as overlap with self:
  out <- as.data.frame(out)
  out2 <- out[- grep(inputlistname, out$listname),]
  
  #add in targetlist name (assumes first list is the input)
  out2$targetlist <- paste(inputlistname)
  
  
  rownames(out2) <- NULL
  
  #return results
  return(out2)
}

##########################################################################################################
#Run Overlap Analysis
##########################################################################################################

#for single list:
#targetlistname <- c("MeCP2 het female MG 24wks")
#genelists <- List
#genomesize <- mm10genome
#MeCP2het24wk <- Overlap_fxn(targetlistname,List,genomesize)



#define targets: list of IDs from List (genelists3)
targets <- names(asdlist)

dfoverlaps <- NULL
for (i in 1:length(targets)) { 
  
  tmp <- Overlap_fxn(targets[i],genelists,genomesize)
  dfoverlaps <- rbind(tmp,dfoverlaps)
  
}

#filter overlaps within experiment:
dfoverlaps <- dfoverlaps[!grepl(".NT", dfoverlaps$listname),]



#adjust pvalue
dfoverlaps$pvalue <- as.numeric(as.character(dfoverlaps$pvalue))
dfoverlaps$FDR <- p.adjust(dfoverlaps$pvalue, method='fdr')

#filter for significant enrichments
SigOverlaps <- dfoverlaps %>% filter(FDR < 0.05)

#Sigcounts$targetlist <- as.factor(Sigcounts$targetlist)
#Sigcounts$listname <- as.factor(Sigcounts$listname)
#Sigcounts <- SigOverlaps %>% dplyr::group_by(targetlist) %>% summarize(count=n())

#filter for overlaps of more than 3 genes:
SigOverlaps$inBinA <- as.numeric(as.character(SigOverlaps$inBinA ))
Sig3plus <- SigOverlaps %>% filter(inBinA > 3)

write.csv(Sig3plus,"SignifcantEnrichments_3plus.csv")
write.csv(dfoverlaps,"allenrichments.csv")


############################################################################################################################################
##########################################################################################################################################
# plots
############################################################################################################################################
############################################################################################################################################

#filter for desired lists:DEGs in ASD, BD or SCZ human brain:
plot <- Sig3plus[grepl(">|<", Sig3plus$listname),]


#get min and max values for heatmaps 
#replace non-significant values with NA

plot$neglog_qvalue <- -log(plot$FDR)

plot$neglog_qvalue[plot$neglog_qvalue == "Inf"] <- 0

plot$neglog_qvalue <- as.numeric(plot$neglog_qvalue)

plot$neglog_qvalue[plot$FDR >=.05] <- NA

#get min and max values for heatmaps 
#sig <- plot2$neglog_qvalue[plot2$FDR <.05]
# 
FCcutoff <- min(plot$neglog_qvalue) -.01
FCcutoff
#[1] 3.193808526
FCcutoff2 <- max(plot$neglog_qvalue,na.rm = T) +.01
FCcutoff2
#[1] 59.441295

#define order:
#order <- plot2 %>% arrange(desc(groups),graphname) %>% dplyr::select(graphname)
#order <- as.character(order$graphname)
#order <- unique(order)

plotPFC <- ggplot(plot, aes(targetlist,listname)) +
  geom_tile(aes(fill = neglog_qvalue),colour = "black") +
  geom_text(aes(label = inBinA))+
  scale_fill_gradient(low = "tan", high = "red",limits = c(FCcutoff,FCcutoff2), guide = "colourbar") +
  scale_y_discrete(name = "Gene List")+ #limits = rev(order)) + #sets orderlimits = nameorder$Term
  xlab(" ") +
  theme(legend.title = element_text(size =15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15,angle = 45, hjust = 1)) +
  labs(fill = "-log(FDR qvalue)")

plotPFC


ggsave("PBMCDEGsoverlapsNDD_DEGsVenn.pdf",width = 8, height = 6, dpi = 300)

#dotplot
#calculate DEG as a percent of each gene list 

percent <- plot[, c('inBinA', 'targetlist')]

a <- vector()
a["Typical.NT>LTA"] <- 278
a["Typical.NT>LPS"] <- 215
a["Typical.NT<LTA"] <- 128
a["Typical.NT<LPS"] <- 122
a["ASD.NT>LTA"] <- 93
a["ASD.NT>LPS"] <- 134
a["ASD.NT<LTA"] <- 86
a["ASD.NT<LPS"] <- 78
a["PDDNOS.NT>LTA"] <- 14
a["PDDNOS.NT>LPS"] <- 47
a["PDDNOS.NT<LTA"] <- 21
a["PDDNOS.NT<LPS"] <- 34

percent$total <- a[percent$targetlist]
percent$percent<- (percent$inBinA / percent$total) * 100
percent$percentround <- round(percent$percent, 1)
plot$percent<- percent$percentround

ggplot(plot, aes(y = targetlist, x = listname)) +
  #facet_grid(~region*DEGdirection, scales = "free")+
  theme(
    panel.grid.major = element_line(colour = "grey"),
    panel.grid.minor = element_line(colour = "grey"),
    legend.title = element_text(size =15),
    legend.text = element_text(size = 15),
    plot.title = element_text(size=15),
    axis.title=element_text(size=15,face="bold"),
    axis.text.y = element_text(size=15),
    axis.text.x = element_text(size=15,angle = 45, hjust = 1))+
  geom_point( alpha=1, aes(size = percent, color=neglog_qvalue)) + 
  xlab(" ") +
  ylab("Overlap List") + 
  scale_color_gradient(low="blue",high="red", limits = c(FCcutoff,FCcutoff2)) 

ggsave("PBMCDEGsoverlapsNDD_DEGsdotpolot.pdf",width = 8, height = 6, dpi = 300)


#############################################################################################################################################################
#Genetic hits 
#############################################################################################################################################################
#filter for desired lists:genetic risk
plot <- Sig3plus[!grepl(">|<", Sig3plus$listname),]



#get min and max values for heatmaps 
#replace non-significant values with NA

plot$neglog_qvalue <- -log(plot$FDR)

plot$neglog_qvalue[plot$neglog_qvalue == "Inf"] <- 0

plot$neglog_qvalue <- as.numeric(plot$neglog_qvalue)

plot$neglog_qvalue[plot$FDR >=.05] <- NA

#get min and max values for heatmaps 
#sig <- plot2$neglog_qvalue[plot2$FDR <.05]
# 
FCcutoff <- min(plot$neglog_qvalue) -.01
FCcutoff
#[1] 3.153412104
FCcutoff2 <- max(plot$neglog_qvalue,na.rm = T) +.01
FCcutoff2
#[1] 14.58076159


plotPFC <- ggplot(plot, aes(targetlist,listname)) +
  geom_tile(aes(fill = neglog_qvalue),colour = "black") +
  geom_text(aes(label = inBinA))+
  scale_fill_gradient(low = "tan", high = "red",limits = c(FCcutoff,FCcutoff2), guide = "colourbar") +
  scale_y_discrete(name = "Gene List")+ #limits = rev(order)) + #sets orderlimits = nameorder$Term
  xlab(" ") +
  theme(legend.title = element_text(size =15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15,angle = 45, hjust = 1)) +
  labs(fill = "-log(FDR qvalue)")

plotPFC


ggsave("PBMCgenetichits_DEGsVenn.pdf",width = 8, height = 3.5, dpi = 300)

#dotplot
#calculate DEG as a percent of target list

percent <- plot[, c('inBinA', 'targetlist')]
percent$total <- a[percent$targetlist]
percent$percent <- (percent$inBinA / percent$total) * 100
percent$percentround <- round(percent$percent, 1)
plot$percent<- percent$percentround

ggplot(plot, aes(y = targetlist, x = listname)) +
  #facet_grid(~region*DEGdirection, scales = "free")+
  theme(
    panel.grid.major = element_line(colour = "grey"),
    panel.grid.minor = element_line(colour = "grey"),
    legend.title = element_text(size =15),
    legend.text = element_text(size = 15),
    plot.title = element_text(size=15),
    axis.title=element_text(size=15,face="bold"),
    axis.text.y = element_text(size=15),
    axis.text.x = element_text(size=15,angle = 45, hjust = 1))+
  geom_point( alpha=1, aes(size = percent, color=neglog_qvalue)) + 
  xlab(" ") +
  ylab("Overlap List") + 
  scale_color_gradient(low="blue",high="red", limits = c(FCcutoff,FCcutoff2)) 

ggsave("PBMCgenetichits_DEGsdotplot.pdf",width = 8, height = 6, dpi = 300)

##########################GO term enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

ensembl <- Sig3plus %>% dplyr::select(intersection_ensembl,listname, targetlist) %>% na.omit
ensembl <- cSplit(ensembl,"intersection_ensembl",",")
ensembl$group <- paste(ensembl$listname,ensembl$targetlist)
ensembl = subset(ensembl, select = -c(1,2))
ensembl <- data.frame(lapply(ensembl, as.character), stringsAsFactors=FALSE)
ensembl_list <- split(ensembl, f = ensembl$group)

enrichmentcc <- compareCluster(geneCluster = ensembl_list,
                               universe      = ASDoverlaplist2$ensembl_gene_id,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "CC",
                               pAdjustMethod = "BH",
                               keyType       = 'ENSEMBL',
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               fun = "enrichGO")

dim(enrichmentcc)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
enrichmentcc <- clusterProfiler::simplify(enrichmentcc, cutoff=0.7, by="p.adjust", select_fun=min)
dim(enrichmentcc)

pdf('humangenes_GO_CCEnrichment.pdf',w=8,h=5)
d <- dotplot(enrichmentcc,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

enrichmentcc <- as.data.frame(enrichmentcc)
write.csv(enrichmentcc,file="humangenes_GO_CCenrichments.csv")

enrichmentbp <- compareCluster(geneCluster = ensembl_list,
                               universe      = ASDoverlaplist2$ensembl_gene_id,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               keyType       = 'ENSEMBL',
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               fun = "enrichGO")

dim(enrichmentbp)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
enrichmentbp <- clusterProfiler::simplify(enrichmentbp, cutoff=0.7, by="p.adjust", select_fun=min)
dim(enrichmentbp)

pdf('humangenes_GO_BPEnrichment.pdf',w=8,h=5)
d <- dotplot(enrichmentbp,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

enrichmentbp <- as.data.frame(enrichmentbp)
write.csv(enrichmentbp,file="humangenes_GO_BP enrichments.csv")

enrichmentmf <- compareCluster(geneCluster = ensembl_list,
                               universe      = ASDoverlaplist2$ensembl_gene_id,
                               OrgDb         = org.Hs.eg.db,
                               ont           = "MF",
                               pAdjustMethod = "BH",
                               keyType       = 'ENSEMBL',
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               fun = "enrichGO")

dim(enrichmentmf)

#simplify: https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
enrichmentmf <- clusterProfiler::simplify(enrichmentmf, cutoff=0.7, by="p.adjust", select_fun=min)
dim(enrichmentmf)

pdf('humangenes_GO_MFEnrichment.pdf',w=8,h=5)
d <- dotplot(enrichmentmf,showCategory = 40)
d + theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

enrichmentmf <- as.data.frame(enrichmentmf)
write.csv(enrichmentmf,file="humangenes_GO_MF enrichments.csv")

##############################

library(enrichR)

dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021","KEGG_2019_Human")
LPS_entrez2 <- subset(mout3, hgnc_symbol %in% LPSDEGs$hgnc_symbol)
LPS_entrez2 <- LPS_entrez2 %>% dplyr::select(LPSVennGroup,hgnc_symbol) %>% na.omit
LPS_entrez2$LPSVennGroup <- gsub("greater"," > ",LPS_entrez2$LPSVennGroup)
TDgLPS <- subset(LPS_entrez2, LPSVennGroup == "Typical.NT > LPS")
enriched <- enrichr(TDgLPS$hgnc_symbol, dbs)
plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
```
