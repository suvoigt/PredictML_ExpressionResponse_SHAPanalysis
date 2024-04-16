#!/bin/bash

# author: Susanne Voigt

help() {
  echo "usage: ./rnaseq_map_count.sh -i <path/to/fastq_ref/folder> -r <path/to/reference.fasta> -g <path/to/reference.gtf>\
  [-t <number_of_threads>][-h <help>]
  - dependencies: STAR, subread 
  - each fastq input file should have a unique name followed by an underscore\
    eg. sample1_R1.fastq.gz and sample1_R2.fastq.gz"
}

# default number of threads = 4 unless specified otherwise by -t option.
threads=4

while getopts "i:r:g:t:h" flag; do
  case "${flag}" in
    i) wdir="${OPTARG}";;
    r) ref="${OPTARG}";;
    g) gtf="${OPTARG}";;
    t) threads="${OPTARG}";;
    h) help; exit 0;;
    *) echo "Invalid option: -$flag" && help && exit 1;;
  esac
done

#-----------
#- mapping -
#-----------

# set working directory
cd "$wdir" || exit


# index reference genome
echo "Indexing reference..."
mkdir -p reference_star_indexed
STAR --runThreadN "$threads" --runMode genomeGenerate \
    --genomeDir reference_star_indexed \
    --genomeFastaFiles "$ref" \
    --sjdbGTFfile "$gtf" --genomeSAindexNbases 12

# list all R1 fastq files 
R1=($(ls -d -1 "$wdir"/* | grep -E '1.fastq(\.gz)?$'))
# list all R2 fastq files 
R2=($(ls -d -1 "$wdir"/* | grep -E '2.fastq(\.gz)?$'))
# get sample names
sample_name=($(ls -1 | grep -E '1.fastq(\.gz)?$' | cut -d'_' -f 1)) # sample names ie. first part of fastq filenames in front of first "_"

echo "Mapping..."
for i in "${!R1[@]}"; do
  unpigz "${R1[$i]}"
  unpigz "${R2[$i]}"
  STAR --genomeDir reference_star_indexed --runThreadN "$threads" \
       --readFilesIn  "${sample_name[$i]}_R1.fastq" "${sample_name[$i]}_R2.fastq" \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts \
       --outFileNamePrefix "${sample_name[$i]}"
  pigz "${sample_name[$i]}_R1.fastq"
  pigz "${sample_name[$i]}_R2.fastq"
done

#-------------------
#- get gene counts -
#-------------------

echo "Getting gene counts..."

# list all gene count output files from STAR
rpg=($(ls -d -1 * | grep -E 'ReadsPerGene.out.tab$'))

# initialize .csv for gene counts
echo -n "" > gene_counts_tmp.csv

# get list of genes
awk 'NR>4 {print $1}' "$(ls *ReadsPerGene.out.tab | head -n 1)" > genes_list_tmp.txt

# iterate over each gene
while read -r gene; do
    echo -n "$gene" >> gene_counts_tmp.csv
    for i in "${rpg[@]}"; do
        # extract the unstranded count for the gene and append to the line
        count=$(awk -v gene="$gene" '$1==gene {print $2}' "$i")
        echo -n ",$count" >> gene_counts_tmp.csv
    done
    echo "" >> gene_counts_tmp.csv
done < genes_list_tmp.txt

# add annotations and header
sort -t, -k1,1 gene_counts_tmp.csv > gene_counts_sorted_tmp.csv
awk '($3=="gene")' "$gtf" | sort -k10,10 | awk 'BEGIN {OFS=","}{print $1,$4,$5,$12,$10}' |
sed 's/"//g' | sed 's/;//g' > gene_annotations_tmp.csv
join -t, -1 5 -2 1 gene_annotations_tmp.csv gene_counts_sorted_tmp.csv > gene_counts_annotations_tmp.csv
echo "gene_id,chrom,start,end,gene_symbol,$(printf "%s," "${sample_name[@]}" | sed 's/,$//')" > header_tmp.csv
cat header_tmp.csv gene_counts_annotations_tmp.csv > gene_counts.csv

rm *_tmp*
