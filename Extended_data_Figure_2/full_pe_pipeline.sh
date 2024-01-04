#!/bin/bash
#PBS -l nodes=1:ppn=3
#PBS -l walltime=48:00:00
#PBS -l pmem=100
#PBS -l feature=rhel7
#PBS -j oe
#PBS -M bur157@psu.edu
#PBS -m bea

#####Catch Error Blocks#############
# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT
#####################################

### FOLDERS #########
DIR=/gpfs/group/adp117/default/bipin/BA_conjugations/2205_RNASeq
RAW=$DIR/00.rawReads
SRC=/gpfs/group/adp117/default/bipin/src

#### SOFTWARES ########
KRAKEN2=/gpfs/group/adp117/default/bipin/src/miniconda3/bin/kraken2

#### DATABASES ##########
KRAKEN_DB=$SRC/databases/kraken2_db
######## INPUT #########
#FILE=$1
FILE=$FILE

##### fastp ########
conda activate rna_seq

FILTERED=$DIR/02.process/a.trim; mkdir -p $FILTERED

FN=$(basename $FILE)

FQ1=$RAW/$FN; SAMPLE=${FN%_R1*}; FQ2=$RAW/${SAMPLE}_R2_001.fastq.gz
FQ1_out=$FILTERED/${SAMPLE}_R1_filtered.fastq.gz; FQ2_out=$FILTERED/${SAMPLE}_R2_filtered.fastq.gz

fastp -i $FQ1 -I $FQ2 -o $FQ1_out -O $FQ2_out -h $FILTERED/${SAMPLE}_fastqc.html  
echo -e "Running fastp on data ${SAMPLE}\n"

conda deactivate

###### Check for contamination #############
conda activate kraken2
echo -e "**Running Kraken and Bracken \n**"


echo -e "Input: $FQ1"
KRAKEN_DIR=$DIR/02.process/d.kraken; mkdir -p $KRAKEN_DIR

$KRAKEN2 --db ${KRAKEN_DB} --gzip-compressed --paired ${FQ1_out} ${FQ2_out} --output ${KRAKEN_DIR}/${SAMPLE}_kraken.out  --report ${KRAKEN_DIR}/${SAMPLE}_kraken.report

conda deactivate

##############################################
conda activate rna_seq 
OUTDIR=$DIR/02.process/b.bam
INDX=$DIR/02.process/c.index/GCF_000196555.1_ASM19655v1_genomic
GTF=$DIR/01.ref/GCF_000196555.1_ASM19655v1_genomic.gtf
OUT=$OUTDIR/$SAMPLE; mkdir -p $OUT
 
bowtie2 -x $INDX -1 $FQ1_out -2 $FQ2_out -S $OUT/${SAMPLE}_aligned.sam

samtools sort $OUT/${SAMPLE}_aligned.sam -o $OUT/${SAMPLE}_aligned.bam

qualimap rnaseq -bam $OUT/${SAMPLE}_aligned.bam -gtf $GTF -pe  -s -outdir $OUT

conda deactivate
# 
# ###############
# 
BAMS=$OUTDIR/bams
find $OUTDIR -iname *.bam -exec cp {} $BAMS/ \;
IN=$BAMS
featureCounts=/gpfs/group/adp117/default/bipin/src/miniconda3/envs/rna_seq/bin/featureCounts
#$featureCounts \
#        -t 'old_locus_tag'\
#        -F 'GTF' \
#	-s 2 \
#	-p --countReadPairs \
#	-g 'CDS' \
#        -a $GTF \
#        -o readCounts.tsv $IN/*.bam

