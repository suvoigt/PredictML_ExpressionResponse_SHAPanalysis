#!/bin/bash

# author: Susanne Voigt

help() {
  echo "usage: ./merge_fastq.sh -i <path/to/info.csv> -d <path/to/folder_with_fastq_files>\

  info.csv should contain a row for each sample with:
  -- col1: sample_name for combined output files
  -- col2: R1_file_1.fastq
  -- col3: R2_file_1.fastq
  -- col4: R1_file_2.fastq
  -- col5: R2_file_2.fastq


  output files: sample_name_R1.fastq.gz and sample_name_R2.fastq.gz"
}


while getopts "i:d:h" flag; do
  case "${flag}" in
    i) csv_file="${OPTARG}";;
    d) wdir="${OPTARG}";;
    h) help; exit 0;;
    *) echo "Invalid option: -$flag" && help && exit 1;;
  esac
done

# set working directory
cd "$wdir"

# read each line from the .csv
while IFS=, read -r combined_name r1_file_1 r2_file_1 r1_file_2 r2_file_2; do
    echo "Processing $combined_name..."

    # concatenate and compress R1 files
    cat "$r1_file_1" "$r1_file_2" | pigz > "${combined_name}_R1.fastq.gz"

    # concatenate and compress R2 files
    cat "$r2_file_1" "$r2_file_2" | pigz > "${combined_name}_R2.fastq.gz"

    # count the lines in each combined file
    r1_lines=$(zcat "${combined_name}_R1.fastq.gz" | wc -l)
    r2_lines=$(zcat "${combined_name}_R2.fastq.gz" | wc -l)

    # check if line counts are equal
    if [ "$r1_lines" -eq "$r2_lines" ]; then
        echo "R1 and R2 combined files for $combined_name have the same number of reads."
        rm "$r1_file_1" "$r1_file_2" "$r2_file_1" "$r2_file_2"
        
    else
        echo "Warning: Mismatch in read counts for $combined_name: R1 has $r1_lines reads, R2 has $r2_lines reads."
    fi
    
    echo "Combined files for $combined_name have been created."
    echo "----------------------------------------------------"
done < "$csv_file"
