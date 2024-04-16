#!/bin/bash

# author: Susanne Voigt

help() {
  echo "usage: ./wgs_trim_map_snps_cns.sh -i <path/to/fastq_ref/folder> -r <path/to/reference.fasta>\
  [-t <number_of_threads>][-h <help>]
  - dependencies: trim_galore, samtools, bwa, gatk, bcftools
  - each fastq input file should have a unique name followed by an underscore\
    eg. sample1_R1.fastq.gz and sample1_R2.fastq.gz"
}

# default number of threads = 4 unless specified otherwise by -t option.
threads=4

while getopts "i:r:t:h" flag; do
  case "${flag}" in
    i) wdir="${OPTARG}";;
    r) ref="${OPTARG}";;
    t) threads="${OPTARG}";;
    h) help; exit 0;;
    *) echo "Invalid option: -$flag" && help && exit 1;;
  esac
done


#----------------------------
#- adapter/quality trimming -
#----------------------------

cd "$wdir"

# make new directory for output
trimmed_dir="$wdir/trimmed"   # directory to trimmed data
mkdir -p "$trimmed_dir"

# list all R1 fastq files 
R1=($(ls -d -1 "$wdir"/* | grep -E '1.fastq(\.gz)?$'))
# list all R2 fastq files 
R2=($(ls -d -1 "$wdir"/* | grep -E '2.fastq(\.gz)?$'))

# quality/adapter trimming using trimgalore
for i in "${!R1[@]}"; do
  trim_galore --paired "${R1[$i]}" "${R2[$i]}" -o "$trimmed_dir" -j "$threads" -q 20 --length 50 --fastqc      # -q min base quality;
done                                                                                                           # --length min read length after trimming


#--------------------------------
#- mapping & quality filtering  -
#--------------------------------

# prepare and indexed reference genome 
awk '{print $1}' "$ref" > ref.fasta                                                  # shorten fasta headers
bwa index ref.fasta                                                                # index reference for bwa
samtools faidx ref.fasta                                                           # index reference for samtools
gatk CreateSequenceDictionary -R ref.fasta -O ref.dict                             # create sequence dictionary for gatk

# Navigate to the directory containing trimmed fastq files


# list all trimmed R1 fastq files
R1_trimmed=($(ls -d -1 "$trimmed_dir"/* | grep -E '1.fq(\.gz)?$'))
# list all trimmed R2 fastq files
R2_trimmed=($(ls -d -1 "$trimmed_dir"/* | grep -E '2.fq(\.gz)?$'))
# get sample names by listing the R1 files and extracting the part before the first "_"
cd "$trimmed_dir"
sample_name=($(ls -1 | grep -E '1.fq(\.gz)?$' | cut -d'_' -f 1))
cd ..

# mapping to reference genome using bwa mem & quality filtering with samtools view
bam_dir="$wdir/bam_vcf"           # directory to bam/vcf data
mkdir -p "$bam_dir"

for i in "${!sample_name[@]}"; do

  bwa mem -M -t "$threads" "ref.fasta" "${R1_trimmed[$i]}" "${R2_trimmed[$i]}" |
  samtools view -Sbh -q 20 -f 0x002 -F 0x104 > "$bam_dir/${sample_name[$i]}_unfiltered.bam"  # -q min mapping quality; only keep proper read pairs(0x002), 
                                                                                                # remove unmapped reads (0x004)  
                                                                                                # and 2ndary alignments (0x100) -> 0x104
  rm "${R1_trimmed[$i]}"
  rm "${R2_trimmed[$i]}"

done


# remove duplicates (with prior sorting by reference position)
for i in "${sample_name[@]}"; do

  # sort bam files by reference position using GATK SortSam
  gatk SortSam \
  -I "$bam_dir/${i}_unfiltered.bam" \
  -O "$bam_dir/${i}_sorted.bam" \
  --SORT_ORDER coordinate
  rm "$bam_dir/${i}_unfiltered.bam"
  
  # mark duplicates and remove them using GATK MarkDuplicates
  gatk MarkDuplicates \
  -I "$bam_dir/${i}_sorted.bam" \
  -O "$bam_dir/${i}_deduped.bam" \
  -M "$bam_dir/${i}_dedup.metrics" \
  --REMOVE_DUPLICATES true
  rm "$bam_dir/${i}_sorted.bam"

done


# get .bam with re-aligned sequences flanking indels and call SNPs (.vcf)
for i in "${sample_name[@]}"; do

  # add read group tags
  gatk AddOrReplaceReadGroups \
    -I "$bam_dir/${i}_deduped.bam" \
    -O "$bam_dir/${i}.bam" \
    -RGID id -RGLB library \
    -RGPL illumina -RGPU unit \
    -RGSM ${i} \
    --CREATE_INDEX true
  rm "$bam_dir/${i}_deduped.bam" 

  # call SNPs and get bam realigned around indels
  gatk HaplotypeCaller  \
   -R ref.fasta \
   -I "$bam_dir/${i}.bam" \
   -O "$bam_dir/${i}.vcf.gz" \
   --native-pair-hmm-threads $threads

done


# mapping statistics
echo -e "\n" > "$bam_dir/mapping_statistics_summary.txt"

for i in "${sample_name[@]}"; do
    # calculate mean depth, mapped reads per chromosome arm
    echo -e "\n$i:\n" >> "$bam_dir/mapping_statistics_summary.txt"
    echo -e "chrom\tchrom_length\tmapped_reads\tmean_depth" >> "$bam_dir/mapping_statistics_summary.txt"
    samtools idxstats "$bam_dir/${i}.bam" |
    awk -v OFS='\t' \
    '($1=="2L"||$1=="2R"||$1=="3L"||$1=="3R"||$1=="4"||$1=="X"||$1=="Y"||$1=="mitochondrion_genome") {if ($2 > 0) {print $1, $2, $3, $3/$2} else {print $1, $2, $3, 0}}' \
    >> "$bam_dir/mapping_statistics_summary.txt"
done



#-----------------------
#- consensus sequence  -
#-----------------------

# get consensus fasta sequence
for i in "${sample_name[@]}"; do
  bcftools view -v snps "$bam_dir/${i}.vcf.gz" -Oz -o "$bam_dir/${i}_snps_only.vcf.gz" # remove indels
  tabix -p vcf "$bam_dir/${i}_snps_only.vcf.gz" # re-index
  bcftools consensus -f ref.fasta "$bam_dir/${i}_snps_only.vcf.gz" -o "${i}.fas" # Generate consensus using the SNP-only VCF
done

