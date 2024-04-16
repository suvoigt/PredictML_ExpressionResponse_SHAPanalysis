library(optparse)
library(edgeR)

# author: Susanne Voigt

# define and parse command-line options:
option_list <- list(
  make_option(
    c('-i', '--input'),
    type = 'character',
    help = 'Path and name of the input file gene annotation(s) followed by gene counts \n
                (.csv, first columns with gene annotations, followed by columns of gene counts per sample)'
  ),
  make_option(
    c('-s', '--samples'),
    type = 'character',
    help = 'Path and name of the file containing sample info (.csv, col1: replicates, col2: group)'
  ),
  make_option(
    c('-o', '--output'),
    type = 'character',
    default = getwd(),
    help = 'Path for output files'
  ),
  make_option(
    c('-a', '--annotation'),
    type = 'integer',
    help = 'Number of columns before the columns with gene counts in input.csv that contain gene annotations\n'
  ),
  make_option(
  c('-c', '--contrast'),
  type = 'character',
  help = 'Comma-separated list of contrasts to perform (eg. "group2 - group1, group3 - group2")'
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


if (!is.null(opt$input)) {
  cat('Input:', opt$input, "\n")

  if (!is.null(opt$samples)) {
    cat('Sample Info:', opt$samples, "\n")

    # read in gene count data
    counts <- read.csv(opt$input)
    # read in sample info
    samples <- read.csv(opt$samples)
    # specify annotaions
    if (!is.null(opt$annotation)) {
        annotations <- counts[, 1:opt$annotation]
        } else {
            cat('Please provide the number of columns containing gene annotations using the --annotation option.\n')
        }
    
    # create DGElist
    dge <- DGEList(counts = counts[, (opt$annotation + 1):ncol(counts)], group = samples[, 2])
    # add annotations to DGElist
    dge$genes <- data.frame(annotations)

    # filter too-lowly expressed genes
    dge_keep <- filterByExpr(dge)
    dge <- dge[dge_keep, keep.lib.sizes=FALSE]

    # calculate normalization factors 
    dge <- calcNormFactors(dge) # default TMM

    # get design matrix - glm analysis
    design_dge <- model.matrix(~0+dge$samples$group)
    colnames(design_dge) <- levels(dge$samples$group)

    # estimate dispersions 
    dge <- estimateGLMCommonDisp(dge,design_dge)
    dge <- estimateGLMTagwiseDisp(dge,design_dge)
    # fit model
    dge_glm <- glmFit(dge,design_dge)


    if (!is.null(opt$contrast)) {
        # parse contrasts
        contrasts <- strsplit(opt$contrast, ',')[[1]]

        # loop through provided contrasts
        for (contrast_name in contrasts) {
        # define contrasts
        contrast <- makeContrasts(
            contrast_name = eval(parse(text = contrast_name)),
            levels = design_dge
        )

        # get differential expression results
        lrt <- glmLRT(dge_glm, contrast = contrast)
        results <- topTags(lrt, n = length(dge$genes$gene_id))

        # write results to .csv
        file_name <- paste0('edger_', gsub("\\s+", "", gsub(" - ", "_vs_", contrast_name)), '.csv')
        write.csv(results, file = file.path(opt$output, file_name))
        } 
        } else {
            cat('Please provide comma-separated list of contrasts (eg. "group2 - group1, group3 - group2") to perform using the --contrast option.\n')
            }

  } else {
    cat('Please provide a file containing information about samples (.csv, col1: replicates, col2: group) using the --samples option.\n')
  }

} else {
  cat('Please provide an input file containing gene annotation(s) followed by gene counts (.csv) using the --input option.\n')
}

