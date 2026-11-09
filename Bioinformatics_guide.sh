## Welcome to Bioinformatics for Buteos! ##

# STEP 1 - Download samples to home directory. Just drop the download.sh file the core facility sent us in the email into your home directory through Filezilla, then run it using this script:

bash download.sh &

#################################################################

# STEP 2 - Run fastqc to look at sequencing quality. You need to run this for each sample. You can either build it here, or make a separate .sh file that lists this command for all individuals, and run it all at once using bash fastq.sh &

fastqc SAMPLENAME.fastq.gz &

#################################################################

# STEP 3 - Adapter removal

#Script for trimming and adapter removal

#Trimns: Trim consecutive Ns from the 5’ and 3’ termini. If quality trimming is also enabled (--trimqualities), then stretches of mixed low-quality bases and/or Ns are trimmed.
#Trimqualities: Trim consecutive stretches of low quality bases (threshold set by --minquality) from the 5’ and 3’ termini. If trimming of Ns is also enabled (--trimns), then stretches of mixed low-quality bases and Ns are trimmed.
#minquality: Set the threshold for trimming low quality bases using --trimqualities and --trimwindows. Default is 2
#collapse: In paired-end mode, merge overlapping mates into a single and recalculate the quality scores. In single-end mode, attempt to identify templates for which the entire sequence is available. In both cases, complete “collapsed” reads are written with a ‘M_’ name prefix, and “collapsed” reads which are trimmed due to quality settings are written with a ‘MT_’ name prefix.
#adapterlist: Read one or more adapter sequences from a table. The first two columns (separated by whitespace) of each line in the file are expected to correspond to values passed to –adapter1 and –adapter2.
#minlength: Reads shorter than this length are discarded following trimming. Defaults to 15Reads shorter than this length are discarded following trimming. Defaults to 15.

# this script uses the default adapters to remove. ALWAYS CHECK TO MAKE SURE YOU DON'T NEED TO DEFINE THE ADAPTER SEQUENCE. If you do, you just add an argument here that specifies a file with the adapter sequence(s) included there. You need to make sure that file is in your working directory if you need it...

# In this format, you need to run the following code for every sample you have. For this and all other times you need to do this, you can either make a separate .sh file that includes the command line for each individual, or figure out how to write a loop.

/programs/adapterremoval_2.1.1/bin/AdapterRemoval --file1 SAMPLENAME_R1.fastq.gz --file2 SAMPLENAME_R2.fastq.gz --trimns --trimqualities --minquality 20 --minlength 25 --collapse --threads 1 --basename SAMPLENAME &

#################################################################

# STEP 4 - Reference genome alignment

# first download reference genome. For this work, we are using a chromosomal level reference for Buteo buteo. When you get to this point, ask either me or bronwyn to put it on the shared directory so you can copy it to your home directory.

#Then you need to run bowtie2-build to create a bowtie index for the reference (you only need to do this once).

bowtie2-build -f Bbuteo.fna Bbuteo

# -f: Reference is in FASTA format. FASTA files usually have extension .fa, .fasta, .mfa, .fna or similar.
#note that the output will be a number of files, they will all have the base name that you specify in NAMEofOUTPUT

#Next, align your reads to the reference genome. You'll need to run the following command line for each sample. I usually do this through independent tmux sessions and fully load a machine.

bowtie2 -p 8 --phred33 --very-sensitive-local \
-x Bbuteo -I 149 -X 900 \
  --rg-id SAMPLENAME --rg SM:SAMPLENAME \
  -1 SAMPLENAME.pair1.truncated -2 SAMPLENAME.pair2.truncated \
  -U SAMPLENAME.collapsed,SAMPLENAME.collapsed.truncated,SAMPLENAME.singleton.truncated \
  2> SAMPLENAME.log | samtools sort -@ 8 -o SAMPLENAME.sorted.bam
samtools index SAMPLENAME.sorted.bam &

# x = NAMEofOUTPUT from step above (base name of your indexed reference)
# p = number of threads (this can be 16 or 24 on the medium machines at BioHPC)
# phred33: refers to the way the quality information is stored in the files - illumina data is --phred33
# --very-sensitive-local is the same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 (these are options that determine the alignment parameters and effort that the program puts into identifying alignments)
# I: The minimum fragment length for valid paired-end alignments. This will depend on the size of the fragments in your library and the length of the reads. If your reads are 150 and are completely overlapping then the min fragement length would be 150.
# X: The maximum fragment length for valid paired-end alignments. This also depends on the size of the fragments in your library and the length of the reads. eg. If your reads are 150bp and are located 600bp apart then they will be considered valid if X=900
#    Note: The larger the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences between -I and -X require that Bowtie 2 scan a larger window to determine if a concordant alignment exists. For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient.
# --rg-id and --rg - these are required for the sam file and adds your SAMPLEID to the read group tag in the sam file. SM means that you are providing the sample name.
# -1 Comma-separated list of files containing first read of mate paired reads (this will be the pair1.truncated file output from AdapterRemoval)
# -2 Comma-separated list of files containing second read of mate paired reads (this will be the pair2.truncated file output from AdapterRemoval)
# -U Comma-separated list of files containing unpaired reads to be aligned (these will be the collapsed, collapsed.truncated and singleton.truncated files output from AdapterRemoval)
# -S File to write SAM alignments to (if you don't specity this your alignments will be written to the standard out file)

#################################################################

# STEP 5 - Run Qualimap to look at mapping statistics. Again, do this for each of your samples. Once they are finished, save all the .pdf files to your local machine in a folder and look through them. There are a number of statistics you need to record in a spreadsheet.

#http://qualimap.conesalab.org/doc_html/index.html

/programs/qualimap_v2.2.1/qualimap bamqc -bam SAMPLENAME_sorted.bam -outfile SAMPLENAME_RT.pdf --java-mem-size=8G &

#################################################################

# STEP 6 - Prep for variant calling

#### Prepare the genome: index it (fai and dict files)
# R:Refernece
# O:Output
java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar CreateSequenceDictionary R=Bbuteo.fna O=Bbuteo.fna.dict

samtools faidx Bbuteo.fna 

# The following will run as a loop. You just need to make sure all of your .sorted.bam files are in the directory that you run the commands in... Look through every command and change the directory and file names where appropriate.

###### mark duplicates - identifies duplicate reads: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
# Metrics File: File name to write duplicate metrics to
# MAX_FILE_HANDLES_FOR_READ_ENDS_MAP: maximum number of file handles to keep open when spilling read ends to a desk. keep this set at 1000
# uses a lot of memory (could probably run more than 8 at a time if I decreased -Xmx to ~10g?)

for RGBAM in *.sorted.bam; do

  #string replacement command
  OUT=${RGBAM/%.sorted.bam/.sorted_mark.bam}
  METRICS=${RGBAM/%.sorted.bam/.metrics.txt}

  #create index file
  echo "java -Xmx15g -jar /programs/picard-tools-2.8.2/picard.jar MarkDuplicates INPUT=$RGBAM OUTPUT=$OUT METRICS_FILE=$METRICS MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000" >> /workdir/bwr46/markDuplicatesCommands.txt

done

parallel -j 8 < /workdir/bwr46/markDuplicatesCommands.txt


# Next you need to validate files

# since running HaplotypeCaller, don't need to realign or fix mates unless there is an error in ValidateSamFile

[ -f /workdir/bwr46/validateCommands.txt ] && rm /workdir/bwr46/validateCommands.txt

cd /workdir/bwr46

for MARKBAM in *mark.bam; do

  OUTFILE=${MARKBAM/%.bam/_validate}

  echo "java -Xmx48g -jar /programs/picard-tools-2.8.2/picard.jar ValidateSamFile I=$MARKBAM OUTPUT=$OUTFILE MODE=SUMMARY" >> /workdir/bwr46/ValidateCommands.txt

done

parallel -j 22 < /workdir/bwr46/ValidateCommands.txt

# concatenate all validate output into one text file to see if there are any errors

cat *validate > summary_validate.txt

#No errors

#### Next index .bam files for HaplotypeCaller
#can run multiple at a time
[ -f /workdir/bwr46/IndexCommands.txt ] && rm /workdir/bwr46/IndexCommands.txt

cd /workdir/bwr46

for MARKBAM in *mark.bam; do

  #create index file
  echo "java -Xmx10g -jar /programs/picard-tools-2.8.2/picard.jar BuildBamIndex I=$MARKBAM" >> /workdir/bwr46/IndexCommands.txt

done

parallel -j 22 < /workdir/bwr46/IndexCommands.txt

#################################################################

# Step 7 - call variants

## CREATE AN ALLSITES VCF THAT INCLUDES FILTERING PERFORMED ON VARIANTS

# Everything should go through PrepForVariantCalling.sh script first...

# Script for variant calling in Samtools
# author: Bryce Robinson. Modified from K Carbeck scripts
# split by chromosome/contig

#### explore a bam file to figure out how to split up over 40 cores
samtools view -H SAMPLENAME.bam_sorted_mark.bam | less

## get genome length
# grep lines that with @SQ
samtools view -H SAMPLENAME.sorted_mark.bam | less | grep "^@SQ" | cut -f 3 | cut -d: -f2 | awk '{sum+=$1}END{print sum}'

## get number of contigs
samtools view -H 10_RT.sorted_mark.bam | grep "^@SQ" | cut -f 3 | cut -d: -f2 | wc -l

# make regions file
samtools view -H 10_RT.sorted_mark.bam | grep "^@SQ" | cut -f 2,3  | awk '
BEGIN {
    cutoff=11640023
    FS=OFS="\t"
    num = 1
}
{
    split($1,ctg,":")
    split($2,ln,":")
    cumulLen += ln[2]
    print ctg[2] , "1" , ln[2] > "regions_" (num<10?"0":"") num ".txt"
    if (cumulLen > cutoff) {
        num++
        cumulLen = 0
    }
}'

# New better way for allsites

bam="SAMPLENAME.sorted_mark.bam"
genome_length=1303682556
n_regions=112
cutoff=$(( (genome_length + n_regions - 1) / n_regions ))

echo "Target size per region: $cutoff bp"




bam="SAMPLENAME.sorted_mark.bam"
n_regions=112

rm -f regions_*.txt

genome_length=$(
    samtools view -H "$bam" |
    awk -F'\t' '
    $1 == "@SQ" {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /^LN:/) {
                split($i, a, ":")
                total += a[2]
            }
        }
    }
    END { print total }
    '
)

cutoff=$(( (genome_length + n_regions - 1) / n_regions ))

echo "Genome length: $genome_length"
echo "Cutoff: $cutoff"

samtools view -H "$bam" |
awk -v cutoff="$cutoff" '
BEGIN {
    FS = OFS = "\t"
    file_num = 1
    file_bases = 0
}

$1 == "@SQ" {
    contig = ""
    contig_length = 0

    for (i = 1; i <= NF; i++) {
        if ($i ~ /^SN:/) {
            contig = substr($i, 4)
        }
        if ($i ~ /^LN:/) {
            contig_length = substr($i, 4) + 0
        }
    }

    start = 1

    while (start <= contig_length) {
        space_left = cutoff - file_bases
        bases_left = contig_length - start + 1

        segment_length = bases_left
        if (segment_length > space_left) {
            segment_length = space_left
        }

        end = start + segment_length - 1
        outfile = sprintf("regions_%03d.txt", file_num)

        print contig, start, end >> outfile
        close(outfile)

        file_bases += segment_length
        start = end + 1

        if (file_bases >= cutoff) {
            file_num++
            file_bases = 0
        }
    }
}
'


# looking at mapping quality
samtools view --exclude-flags SECONDARY,SUPPLEMENTARY SAMPLE_sorted_mark.bam | cut -f5 | sort -nr | uniq -c | less


# http://www.htslib.org/doc/bcftools.html#mpileup
# Generate VCF or BCF containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files. This is based on the original samtools mpileup command (with the -v or -g options) producing genotype likelihoods in VCF or BCF format, but not the textual pileup output. The mpileup command was transferred to bcftools in order to avoid errors resulting from use of incompatible versions of samtools and bcftools when using in the mpileup+bcftools call pipeline.
# Individuals are identified from the SM tags in the @RG header lines. Multiple individuals can be pooled in one alignment file, also one individual can be separated into multiple files. If sample identifiers are absent, each input file is regarded as one sample.

## mpileup
  # -O: output type. u = bcf
  # -C: Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50
  # -a: annotate. list of common INFO tags for vcf output
  # -f: Reference genome
## call
  # -multiallelic-caller: alternative model for multialleleic and rare-variant calling designed to overcome known limitations in -c calling mode (previous version)
  # -variants-only: output variant sites only
  # -O: output type. v=VCF
  # -o: output file name
  # -b = List of input alignment files, one file per line (samplename_sorted_RGmark.bam)


mkdir vcf2


## Need to remove variants only flag plus any -v

genome=Bbuteo.fna
bamfiles=*.sorted_mark.bam
for regions in regions_*.txt; do
   mpileupVCF=${regions%.txt}.raw.mpileup.vcf
    echo $mpileupVCF
    bcftools mpileup -a "FORMAT/AD,FORMAT/DP,INFO/AD,FORMAT/DV,FORMAT/DP4,FORMAT/DPR,INFO/DPR" --threads 2 -d 1000 --min-MQ 30 --regions-file $regions -Ou -f $genome $bamfiles | bcftools call --multiallelic-caller --threads 2 -m -Ov -o vcf2/$mpileupVCF 2>&1 | tee vcf2/${regions%.txt}.mpileup.$(date +%Y%m%d-%Hh%Mm%Ss).log &
    sleep 30
done


#################################################################################################
#### bcftools concatenate ####
# you can't just cat them together because each file has a header section. The bcftools command will handle that.

# first, each vcf file must be sorted prior to calling bcftools concat
mkdir -p vcf2/sorted

max_jobs=12

for regions in regions_*.txt; do

    prefix=$(basename "${regions}" .txt)

    mpileupVCF="vcf2/${prefix}.raw.mpileup.vcf"
    sortedVCF="vcf2/sorted/${prefix}.sorted.vcf.gz"

    if [[ ! -s "${mpileupVCF}" ]]; then
        echo "ERROR: missing or empty input: ${mpileupVCF}" >&2
        continue
    fi

    echo "Sorting ${mpileupVCF} -> ${sortedVCF}"

    bcftools sort \
        "${mpileupVCF}" \
        -Oz \
        -o "${sortedVCF}" &

    while [[ $(jobs -rp | wc -l) -ge ${max_jobs} ]]; do
        sleep 2
    done

done

wait

# then index

find vcf2/sorted \
    -name "regions_*.sorted.vcf.gz" \
    -type f |
    wc -l

# then concat

find vcf2/sorted \
    -name "regions_*.sorted.vcf.gz" \
    -type f |
    sort -V > concat_vcf_list.txt

bcftools concat \
    -f concat_vcf_list.txt \
    -Oz \
    -o vcf2/raw_ALLSITES.vcf.gz

# index again??

bcftools index -t vcf2/raw_ALLSITES.vcf.gz



# valadate vcf to make sure everything is ok

vcf-validator raw_ALLSITES.vcf.gz

#15,767,149 variants


# filter allsites

#!/bin/bash
set -euo pipefail

vcftools \
    --gzvcf raw_ALLSITES.vcf.gz \
    --min-meanDP 5 \
    --max-meanDP 80 \
    --max-missing 0.8 \
    --minQ 35 \
    --recode --stdout |
bgzip -c > filtered_ALLSITES.vcf.gz

tabix -f -p vcf filtered_ALLSITES.vcf.gz

# Creat variants only

bcftools view \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --min-af 0.05:minor \
    -Oz \
    -o filtered_variants_MAF05.vcf.gz \
    filtered_ALLSITES.vcf.gz

bcftools index -f -t filtered_variants_MAF05.vcf.gz

# Record number of variants - 7,272,375

bcftools view -H filtered_variants_MAF05.vcf.gz | wc -l

# list indivduals to make sure everyone is included

bcftools query -l filtered_variants_MAF05.vcf.gz

###################################################
## ADMIXTURE ##
###################################################

#Grasshopper Sparrow Project
#Jen Walsh, March 2023
#Analyses - ADMIXTURE

#Remove missing data and prune for linkage

vcftools --gzvcf filtered_variants_MAF05.vcf.gz --max-missing 1 --recode --stdout | bgzip > filtered_variants_MAF05.vcf.gz

#Prune for linkage
bash ldPruning.sh filtered_variants_MAF05.vcf.gz

#Record number of snps retained:

#Prepare input files for ADMIXTURE
vcftools --gzvcf filtered_variants_MAF05.LDpruned.vcf.gz --plink --out RTHA.noN.LDprunedPlinkformat
plink --file RTHA.noN.LDprunedPlinkformat --make-bed

# 200 iterations for k 1-4 (which is the number of sites/populations sampled)
for K in 1 2 3 4; \
	do admixture --cv=5 plink.bed $K | tee RTHA_K${K}.out; done

# -C 100
#^flag to set the number of iterations to 100

echo "CV Error:"
grep -h CV RTHA_K*.out

#CV Error:
#CV error (K=1): 
#CV error (K=2): 
#CV error (K=3): 
#CV error (K=4): 
#CV error (K=5): 
#CV error (K=6): 

# Then plot in R

