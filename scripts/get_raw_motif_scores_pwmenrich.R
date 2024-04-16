library(optparse)
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)


# author: Susanne Voigt


# define and parse command-line options:
option_list <- list(
  make_option(
    c('-i', '--input'),
    type = 'character',
    help = 'Path and name of the input file containing dataframe (col1: gene, col2: path/to/fasta)'
  ),
  make_option(
    c('-o', '--output'),
    type = 'character',
    default = 'output_raw_scores_pwmenrich.csv',
    help = 'Path and name of the output file (default: "output_raw_scores_pwmenrich.csv")'
  ),
  make_option(
    c('-n', '--n_group1'),
    type = 'integer',
    help = 'Number of sequences of group 1.\n
                Multiple-sequence-alignment fasta files must consist first of the sequences of group 1, followed by those of group 2.\n 
                If not provided, motif scores are calculated across all sequences of fasta.'
  ),
  make_option(
    c('-c', '--cores'),
    type = 'integer',
    default = 1,
    help = 'Number of cores to utilize (default: 1)'
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)



# functions:
# calculate raw motif scores for 1-group fasta files
motifScores_1 <- function(fas, gene) {
  # read in fasta
  seq <- readDNAStringSet(fas) 
  # calculate motif scores for each group of sequences
  res <- motifEnrichment(seq, MotifDb.Dmel_inclDsp1)
   # create dataframe
  df <- data.frame(motif=groupReport(res1)$id, 
                    score1=groupReport(res1)$raw.score) # extract motif IDs and scores
  # add gene-specific column names
  colnames(df)[2] <- paste0(gene, '_score')
  # return scores
  return(df)
}

# calculate raw motif scores for 2-groups fasta files 
motifScores_2 <- function(fas, gene, n_group1) {
  # read in fasta
  seq <- readDNAStringSet(fas) 
  # split sequences by group
  seq1 <- seq[1:n_group1]
  seq2 <- seq[(n_group1+1):length(seq)]
  # calculate motif scores for each group of sequences
  res1 <- motifEnrichment(seq1, MotifDb.Dmel_inclDsp1)
  res2 <- motifEnrichment(seq2, MotifDb.Dmel_inclDsp1)
  # create dataframe for each group and sort by motif id
  df1 <- data.frame(motif=groupReport(res1)$id, 
                    score1=groupReport(res1)$raw.score) # extract motif IDs and scores
  df2 <- data.frame(motif=groupReport(res2)$id, 
                     score2=groupReport(res2)$raw.score) # extract motif IDs and scores
  # combine the scores of both groups into one dataframe
  df <- merge(df1, df2, by='motif', all=TRUE)
  # add gene-specific column names
  colnames(df)[2:3] <- c(paste0(gene, '_score1'), paste0(gene, '_score2'))
  # return scores
  return(df)
}


# main:
if (!is.null(opt$input)) {
  cat('Input:', opt$input, "\n")
  
  # parallel excution 
  registerCoresPWMEnrich(opt$cores)

  # load motif data
  data(MotifDb.Dmel)
  # remove motifs in MotifDb.Dmel with duplicated names
  MotifDb.Dmel <- MotifDb.Dmel[-(which(duplicated(names(MotifDb.Dmel), fromLast = TRUE)))]
  # add Dsp1 motif to MotifDb.Dmel (retrieved from Jaspar https://jaspar.genereg.net/matrix/UN0108.1/):
  motif.Dsp1 <- readMotifs('data/meta/Dsp1.jaspar',remove.acc=TRUE)
  # convert count matrices into PWMs
  genomic.acgt <- getBackgroundFrequencies('dm3')
  pwms.Dsp1 <- toPWM(motif.Dsp1, prior=genomic.acgt)
  #append pwms.Dsp1 to MotifDb.Dmel
  MotifDb.Dmel_inclDsp1 <- append(MotifDb.Dmel, pwms.Dsp1) 

  # read in dataframe with input specifications
  df <- read.csv(opt$input)
  # extract lists of gene names and paths to fasta files from dataframe
  gene <- df$gene_symbol
  fas <- df$fasta

  if (is.null(opt$n_group1)) {
    # based on first fasta - set up dataframe for results with first columns providing info about TFs and motifs
    seq1 <- readDNAStringSet(fas[1]) 
    res1 <- motifEnrichment(seq1, MotifDb.Dmel_inclDsp1) # calculate motif scores
    scores1 <- data.frame(TF=groupReport(res1)$target, 
                          motif=groupReport(res1)$id, 
                          score1=groupReport(res1)$raw.score) # extract motif IDs and scores into df
    colnames(scores1)[3] <- paste0(gene[1], '_score1') # rename columns with scores for specific gene region

    # run motifScores function on all other gene regions
    scores <- mapply(motifScores_1,fas[2:length(fas)], gene[2:length(fas)], SIMPLIFY=FALSE)
    # combine list of dataframes in one dataframe
    scores_df <- Reduce(function(x, y) merge(x, y, by='motif', all=TRUE), scores)

    # merge dataframes for output
    output_df <- merge(scores1, scores_df, by='motif', all=TRUE)
    # create output file
    write.csv(output_df, opt$output)

  } else {
     # based on first fasta - set up dataframe for results with first columns providing info about TFs and motifs
    seq1 <- readDNAStringSet(fas[1]) # read in fasta
    seq1_1 <- seq1[1:opt$n_group1] # split sequences by group 1 and 2 - group 1 sequences
    seq1_2 <- seq1[(opt$n_group1+1):length(seq1)] # split sequences by group 1 and 2 - group 2 sequences
    res1_1 <- motifEnrichment(seq1_1, MotifDb.Dmel_inclDsp1) # calculate motif scores - group 1
    res1_2 <- motifEnrichment(seq1_2, MotifDb.Dmel_inclDsp1) # calculate motif scores - group 2
    scores1_1 <- data.frame(TF=groupReport(res1_1)$target, 
                            motif=groupReport(res1_1)$id, 
                            score1=groupReport(res1_1)$raw.score) # extract motif IDs and scores into df - group 1
    scores1_2 <- data.frame(motif=groupReport(res1_2)$id, 
                            score2=groupReport(res1_2)$raw.score) # extract motif IDs and scores into df - - group 2
    scores1 <-  merge(scores1_1, scores1_2, by='motif', all=TRUE) # merge dfs of group1 and group2 scores 
    colnames(scores1)[3:4] <- c(paste0(gene[1], '_score1'), paste0(gene[1], '_score2')) # rename columns with scores for specific gene region

    # run motifScores function on all other gene regions
    scores <- mapply(motifScores_2,fas[2:length(fas)], gene[2:length(fas)], opt$n_group1, SIMPLIFY=FALSE)
    # combine list of dataframes in one dataframe
    scores_df <- Reduce(function(x, y) merge(x, y, by='motif', all=TRUE), scores)

    # merge dataframes for output
    output_df <- merge(scores1, scores_df, by='motif', all=TRUE) 
    # create output file
    write.csv(output_df, opt$output)
  }

} else {
  cat('Please provide an input file using the --input option.\n')
}