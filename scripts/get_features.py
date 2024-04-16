import os
import sys
import argparse
import pandas as pd
import numpy as np


# author: Susanne Voigt


# functions 
def check_file_path(file_path):
    '''
    Checks if the specified file path exists.

    Args:
      file_path (str): The path to the file.

    Returns:
      str: The validated file path if it exists.

    Raises:
      argparse.ArgumentTypeError: If the file path does not exist.
    '''
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"{file_path} does not exist.")
    return file_path


def clean_up_gene_names(df):
    '''
    Cleans up gene symbol names in a DataFrame.

    Args:
      df (DataFrame): DataFrame containing a column 'gene_symbol' with gene names.

    Returns:
      DataFrame: The input DataFrame with cleaned up 'gene_symbol' column.
    '''
    # remove genes on chrom 4, rename gene 'nan' to full name 'nanchung', and replace ':' in gene symbols with '_'
    df = df.loc[df['chrom'] != '4', :].reset_index(drop=True)
    df['gene_symbol'] = df['gene_symbol'].fillna('nanchung')
    df['gene_symbol'] = df['gene_symbol'].str.replace(':', '-')

    return df


def get_motif_score_features(motif_data, group, motif_info, annotation):
    '''
    Processes motif data to extract, clean, and organize motif score features for a specified group and 
    sets a detailed index based on genomic location and gene symbol from an external `annotation` DataFrame.

    Args:
      motif_data (DataFrame): A DataFrame containing motif scores and identifiers.
      group (int): The specific group number to extract motif score features for.
      motif_info (DataFrame): DataFrame containing motif IDs and associated information.
      annotation (DataFrame): DataFrame containing annotations of considered genes.

    Returns:
      DataFrame: A DataFrame with motif score features for the specified group, sorted by index,
                 with a custom index based on genomic location, gene symbol, and group number.
    '''
    # clean motif data
    df = pd.merge(motif_info, motif_data, left_on='motif_ID', right_on='motif', how='inner')
    df = df.drop(df.columns[0:4], axis=1).drop(df.columns[5:7], axis=1)
    
    # extract group-specific scores and rename columns
    df = df.loc[:, ['TF_motif'] + [col for col in df.columns if f'_score{group}' in col]]
    df['TF_motif'] = df['TF_motif'] + '_score'
    df = df.rename(columns={'TF_motif': 'pos_gene_group'})
    
    # transpose and set custom index
    df = df.set_index('pos_gene_group').T.sort_index()
    df_genes_anno = annotation[annotation['gene_symbol'].isin(df.index.str.replace(f'_score{group}', ''))].sort_values(by='gene_symbol')
    df.index = df_genes_anno['chrom'] + '_' + df_genes_anno['start'].astype(str) + '_' + \
               df_genes_anno['end'].astype(str) + '_' + df_genes_anno['gene_symbol'] + '_' + str(group)
    
    return df


def get_expr_features(expr_data, group, annotation):
    '''
    Transforms expression data (edgeR output file) into a DataFrame of gene expression features. 

    Args:
        expr_data (DataFrame): Contains 'gene_symbol' and 'logFC'.
        group (int): The specific group number to extract motif score features for.
        annotation (DataFrame): DataFrame containing annotations of considered genes.

    Returns:
        DataFrame: Each column represents a gene's logFC, prepared for analysis.
   '''
    # extract fold-changes for each gene
    df = expr_data.loc[:, ['gene_symbol', 'logFC']]
    df['gene_symbol'] += '_expr_logFC'

    # transpose, expand and set custom index
    df = df.rename(columns={'gene_symbol': 'pos_gene_group'})
    df = df.set_index('pos_gene_group').sort_index().T
    df = df.loc[np.repeat(df.index, len(annotation))].reset_index(drop=True)
    df.index = annotation['chrom'] + '_' + annotation['start'].astype(str) + '_' + \
               annotation['end'].astype(str) + '_' + annotation['gene_symbol'] + '_' + str(group)

    return df


def get_target(expr_data, group, annotation):
    '''
    Transforms expression data (edgeR output file) into a DataFrame of the binary target vector
    'TSE', ie. significant overexpression at lower temperature with '1' for significant overexpression 
    at lower temperature (P<0.05 & logFC>0), and '0' for none. 

    Args:
        expr_data (DataFrame): DataFrame containing expression data.
        group (int): The specific group number to extract motif score features for.
        annotation (DataFrame): DataFrame containing annotations of considered genes.

    Returns:
        df (DataFrame): One column 'TSE' of target vector, ie. significant overexpression at lower temperature.
   '''
    # extract expression of target genes, and encode into binary
    expr_target_genes = expr_data[expr_data['gene_symbol'].isin(annotation['gene_symbol'])].sort_values(by='gene_symbol').reset_index()
    expr_target_genes['TSE'] = expr_target_genes.apply(lambda row: 1 if row['PValue']<0.05 and row['logFC']>0 else 0, axis=1)
    df = expr_target_genes[['TSE']].copy()

    # set custom index
    df_genes_anno = annotation[annotation['gene_symbol'].isin(expr_target_genes['gene_symbol'])].sort_values(by='gene_symbol').reset_index()
    df['pos_gene_group'] = df_genes_anno['chrom'] + '_' + df_genes_anno['start'].astype(str) + '_' + \
                           df_genes_anno['end'].astype(str) + '_' + df_genes_anno['gene_symbol'] + '_' + str(group)
    df.set_index('pos_gene_group', inplace=True) 

    return df


def get_baseline(expr_data, group, annotation):
    '''
    Transforms expression data (edgeR output file) into a DataFrame containing the mean gene activity 
    (log Counts Per Million, logCPM) across all samples for use as a baseline predictor.

    Args:
        expr_data (DataFrame): DataFrame containing expression data.
        group (int): The specific group number to extract motif score features for.
        annotation (DataFrame): DataFrame containing annotations of considered genes.

    Returns:
        df (DataFrame): One column with mean gene activity (ie., mean logCPM across all samples)
    '''
    # extract gene activity
    expr_target_genes = expr_data[expr_data['gene_symbol'].isin(annotation['gene_symbol'])].sort_values(by='gene_symbol').reset_index()
    df = expr_target_genes[['logCPM']].copy()

    # set custom index
    df_genes_anno = annotation[annotation['gene_symbol'].isin(expr_target_genes['gene_symbol'])].sort_values(by='gene_symbol').reset_index()
    df['pos_gene_group'] = df_genes_anno['chrom'] + '_' + df_genes_anno['start'].astype(str) + '_' + \
                           df_genes_anno['end'].astype(str) + '_' + df_genes_anno['gene_symbol'] + '_' + str(group)
    df.set_index('pos_gene_group', inplace=True) 


    return df


if __name__ == "__main__":

    # define and parse command-line options:
    parser = argparse.ArgumentParser(description='A script to build feature dataframe of TFBS motif scores  of one or more groups/lines/populations \
                                                  (with each row corresponding to a specific gene (region). \
                                                  If expression data is provided, the dataframe will also include expression features (logFCs), \
                                                  a binary target vector (ie. significant overexpression in positive direction (1) or none (0), \
                                                  as well as as a baseline predictor dataset consisting of \
                                                  gene activity (ie, mean logCPM values per gene) only feature-target dataframe.')

    parser.add_argument('--motif_scores', type=check_file_path, required=True, help='Path to get_raw_motif_scores_pwmenrich.R output .csv with raw motif scores.')
    parser.add_argument('--motif_info', type=check_file_path, required=True, help='Path to .csv with information about non-redudandant motifs (col1:TF, col2: motif_ID, col3: consensus, col4: data_source, col5: TF_motif).')
    parser.add_argument('--annotation', type=check_file_path, required=True, help='Path to .csv with annotations of genes (col1: gene_id, col2: chrom, col3: start, col4: end, col5: gene_name).')
    parser.add_argument('--n_groups', type=int, default=1, required=True,help='Number of separate groups/lines/populations.')
    parser.add_argument('--expression', type=check_file_path, nargs='*', help='Optional: Path to to edgeR output .csv file(s) with genome-wide expression results. You can specify multiple files (eg. --expression file_group1.csv file_group2.csv) corresponding to n_groups option.')
    parser.add_argument('--output', default='features.csv', help='Path to output file. Defaults to "features.csv" in the current working directory.')
    parser.add_argument('--output2', default='gene_activity.csv', help='Optional: Path to output file of baseline predictor. Defaults to "gene_activity.csv" in the current working directory.')

    args = parser.parse_args()

    # validate the number of expression files against n_groups if expression files are provided
    if args.expression and len(args.expression) != args.n_groups:
        print(f'Error: The number of expression files provided ({len(args.expression)}) does not match the expected number of groups ({args.n_groups}).', file=sys.stderr)
        sys.exit(1)

    # main
    motif_scores = pd.read_csv(args.motif_scores)
    motif_info = pd.read_csv(args.motif_info)
    annotation = clean_up_gene_names(pd.read_csv(args.annotation))

    motif = [get_motif_score_features(motif_scores, n, motif_info, annotation) for n in range(1, args.n_groups+1)]
    motif_features = pd.concat(motif,ignore_index=False)

    if args.expression:
        # expression data provided
        print('Feature dataframe including motif scores and expression data, plus target vector will be built, as well as a baseline predictor dataset (ie., mean gene activity plus target vector).')
        
        expr = []
        target = []
        base = []
        for i, file_path in enumerate(args.expression, start=1):
            expr_data = pd.read_csv(file_path)
            expr_data_clean = clean_up_gene_names(expr_data)
            expr.append(get_expr_features(expr_data_clean, i, annotation))
            target.append(get_target(expr_data_clean, i, annotation))
            base.append(get_baseline(expr_data_clean, i, annotation))
        
        expr_features = pd.concat(expr,ignore_index=False)
        target_vector = pd.concat(target,ignore_index=False)
        baseline = pd.concat(base,ignore_index=False)

        expr_features_target = pd.merge(expr_features, target_vector, left_index=True, right_index=True, how='inner')
        features_target = pd.merge(motif_features, expr_features_target, left_index=True, right_index=True, how='inner')
        features_target.to_csv(args.output)

        baseline_target = pd.merge(baseline, target_vector, left_index=True, right_index=True, how='inner').loc[features_target.index]
        baseline_target.to_csv(args.output2)


    else:
        # no expression data provided
        print('No expression files provided. Feature dataframe including motif scores only will be built.')

        motif_features.to_csv(args.output)




