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
