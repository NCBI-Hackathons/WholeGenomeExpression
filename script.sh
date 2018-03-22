#!/bin/bash
# Usage: sh script.sh SRR531311 SRR531315 ./RefDir 8 ./outDir > log &

if [ "$#" -ne 5 ]; then
       echo "Not Enough Arguments. Please enter sh script.sh [SRA ID #1] [SRA ID #2] [Reference Genome Dir] [CPUs] [Output Directory]"
       exit 1
fi

SRA1ID="$1"
SRA2ID="$2"
REFDIR=$(readlink -f "$3")
CORES="$4"
OUTDIR=$(readlink -f "$5")

#printf "Getting host genome and create BLAST databases\n"
#date
#mkdir $hostGen
#wget -P $hostGen ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
#gunzip $hostGen/GRCh37_latest_genomic.fna.gz
#makeblastdb -in $hostGen/GRCh37_latest_genomic.fna -dbtype nucl -out $hostGen/hg19

printf "Getting Reference genome and create BLAST databases\n"
date
mkdir $REFDIR
wget -P $REFDIR `esearch -db assembly -query 'Mus musculus [orgn] AND latest [SB]' | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_genomic.fna.gz"}'` 
gunzip $REFDIR/*.fna.gz
makeblastdb -in $REFDIR/*_genomic.fna -dbtype nucl -out $REFDIR/mouseGenome
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
wget -P $REFDIR ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz
gunzip $REFDIR/GCF_000001635.26_GRCm38.p6_genomic.gff.gz

mkdir $OUTDIR
printf "Creating SAM alignment files with magicblast\n"
date
magicblast -sra $SRA1ID -db $REFDIR/mouseGenome -out $OUTDIR/$SRA1ID.mm.sam -num_threads $CORES -no_unaligned
magicblast -sra $SRA2ID -db $REFDIR/mouseGenome -out $OUTDIR/$SRA2ID.mm.sam -num_threads $CORES -no_unaligned

printf "Convert SAM to BAM\n"
date
#samtools view -bS $OUTDIR/$SRA1ID.mm.sam | samtools sort -o $OUTDIR/$SRA1ID.mm.bam
#samtools view -bS $OUTDIR/$SRA2ID.mm.sam | samtools sort -o $OUTDIR/$SRA2ID.mm.bam
samtools view -bS $OUTDIR/$SRA1ID.mm.sam | samtools sort - $OUTDIR/$SRA1ID.mm
samtools view -bS $OUTDIR/$SRA2ID.mm.sam | samtools sort - $OUTDIR/$SRA2ID.mm

printf "Counting\n"
date
featureCounts -t exon -g gene -a $REFDIR/GCF_000001635.26_GRCm38.p6_genomic.gff -o counts_gene1.txt $OUTDIR/$SRA1ID.mm.bam
featureCounts -t exon -g gene -a $REFDIR/GCF_000001635.26_GRCm38.p6_genomic.gff -o counts_gene2.txt $OUTDIR/$SRA2ID.mm.bam

printf "Done\n"
date
