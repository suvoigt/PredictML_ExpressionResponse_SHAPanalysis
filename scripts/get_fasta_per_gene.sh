#!/bin/bash

# author: Susanne Voigt

help() {
  echo "usage: ./get_fasta_per_gene.sh -i <path/to/input1.fas> [-i <path/to/input2.fas> ...] -r <path/to/regions.bed> [-h <help>]
  - dependencies: bedtools, mafft
  - regions.bed defines the segments to extract from input FASTA files
    (tab-separated columns with col1: chrom, col2: start, col3: end, col4: name)"
}

input_files=() # Array to hold input FASTA files

while getopts "i:r:h" flag; do
  case "${flag}" in
    i) input_files+=("${OPTARG}");;
    r) bed="${OPTARG}";;
    h) help; exit 0;;
    *) echo "Invalid option: -$flag" && help && exit 1;;
  esac
done

if [ ${#input_files[@]} -eq 0 ] || [ -z $bed ]; then
  echo "Error: input .fasta files and regions.bed are required."
  help
  exit 1
fi

# make new directory for output
output_dir=$(dirname ${input_files[0]})/fasta_per_gene
mkdir -p $output_dir

# Process each region defined in .bed
while IFS= read -r line || [[ -n $line ]]; do
  read -ra arr <<<$line
  chrom=${arr[0]}
  start=${arr[1]}
  end=${arr[2]}
  gene=${arr[3]}

  # file to compile sequences from all inputs for this region
  compiled_seqs="${output_dir}/${gene}_${chrom}_${start}_${end}_compiled.fas"

  # ensure the compiled file is empty or does not exist
  > $compiled_seqs

  # extract sequences from each input for the current region and concatenate them
  for input in ${input_files[@]}; do
    # extract region and append to the compiled file
    echo -e "${chrom}\t${start}\t${end}" | bedtools getfasta -fi $input -bed - -fo - -name |
    sed "s/>/>$(basename "$input" .fas)/" >> $compiled_seqs
  done

  # align sequences using MUSCLE and output to a new file
  aligned_seqs="${compiled_seqs/_compiled.fas/.fas}"
  mafft "$compiled_seqs" | awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {if (seq) print seq}' |
  awk '{if(/^>/) {print} else {print toupper(gensub(/-/, "N", "g", $0))}}' > "$aligned_seqs"
  rm "$compiled_seqs"


done < $bed

