# Kenneth B. Hoehn, Mamie Wang
# 11/22/2021
# Script to cluster sequences by heavy/ heavy+light chain by multiple methods
# using Scoper
#
# Input:
# AIRR-tsv formatted table of heavy and light chain genes
# required columns: sequence_id, sequence_alignment, cdr3, v_call, j_call, junction,
# cell_id, locus; germline_alignment required for spectral vj
#
# Output:
# AIRR-tsv formatted table with clone ids in the clone_id column
#
# Run:
# Rscript scoperClones.R <input.tsv> <output.tsv> <hierarchical|spectral> <H|HL> <novj|vj> <nproc> <thresh>
# 
# Options:
# hierarchical: Use single-linkage hierarhical clustering with automated threshold detection
# spectral: Use spectral clustering
# H: Use only heavy chain
# HL: Split clones by light chain
# novj: Don't use shared mutation information
# vj: Use shared mutation information (only possible with "spectral")
# nproc: Number of cores for parallel processing
# thresh: optional argument to manually specify cutoff threshold for hierarchical clustering
#
# Notes:
#    - Cell with multiple heavy chains will be removed
#    - clone_id of cells with only light chains will be NA
#
# References
# vj/novj: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007977
# hierarhical clustering: https://www.jimmunol.org/content/198/6/2489
# spectral clustering: https://academic.oup.com/bioinformatics/article/34/13/i341/5045726

library(alakazam)
library(shazam)
library(dplyr)
library(scoper)

input = "combined.tsv"
output = "output.tsv"
clustering = "hierarchical"
chain = "HL"
vj = "novj"
nproc = 1

args = commandArgs(trailingOnly=TRUE)

argsLen <- length(args);
if (argsLen > 7) stop('error: too many arguments.');
if (argsLen < 6) stop('error: missing arguments')

input = args[1]
output = args[2]
clustering = args[3]
chain = args[4]
vj = args[5]
nproc = as.numeric(args[6])
if (argsLen == 7) {
    threshold = as.numeric(args[7])
}

data = readChangeoDb(input)

# Check options
if(chain == "HL"){
	split_light = TRUE
}else if(chain == "H"){
	split_light = FALSE
}else{
	stop("Invalid chain option (must be H or HL")
}
if(!vj %in% c("vj","novj")){
	stop("Invalid vj specification (must be novj or vj)")
}
if(vj == "vj" && clustering != "spectral"){
	stop("vj option only available with spectral clustering")
}

# Remove cells with multiple heavy chains
multi_heavy = table(filter(data, locus=="IGH")$cell_id)
multi_heavy_cells = names(multi_heavy)[multi_heavy > 1]
data = filter(data, !cell_id %in% multi_heavy_cells)
if(length(multi_heavy_cells) > 0){
	print(paste("Removed", length(multi_heavy_cells), "cells with multiple heavy chains"))
}

if(clustering == "hierarchical"){    
    # find clonal thresholds    
    # https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/
    # uses cdr3 instead of junction to be compatible with cdr3=TRUE in hierarchichalClones
    # can change sequenceColumn to "junction" and cdr3=FALSE for hierarchicalClones to 
    # use the full junction region instead
    # also, threshold is only calculated using the heavy chain
    dist_ham <- distToNearest(filter(data, locus=="IGH"),     
        sequenceColumn="cdr3", 
        vCallColumn="v_call",
        jCallColumn="j_call",
        model="ham", 
        normalize="len",
        nproc=nproc)
    if (argsLen == 6) {
        output <- findThreshold(dist_ham$dist_nearest, method="density")
        threshold <- output@threshold
    }
    

    # uncomment if threshold plot desired
    #    pdf("Threshold.pdf",height=6,width=6)
    #    plotDensityThreshold(output)
    #    dev.off()

    # stop if threshold detection fails
    if(is.na(threshold)){
        stop("Threshold detection failed, impossible to continue with hierarchical clustering")
    }

    # perform hierarchical clustering with specified options
    # https://scoper.readthedocs.io/en/stable/topics/hierarchicalClones/    
    results = hierarchicalClones(data,
        threshold=threshold,
        cdr3=TRUE,
        only_heavy=!split_light,
        cell_id="cell_id",
        locus="locus",
        split_light=split_light,
        nproc=nproc
        )
    results = as.data.frame(results)

}else if(clustering == "spectral"){
	
	# perform spectral clustering with specified options
	# https://scoper.readthedocs.io/en/stable/topics/spectralClones/
    results = spectralClones(data,
        cdr3=TRUE,        
        only_heavy=!split_light,
        cell_id="cell_id",
        locus="locus",
        method=vj,
        split_light=split_light,
        nproc=nproc
        )
    results = as.data.frame(results)

}else{    
	stop("Invalid clustering specification (must be spectral or hierarchical)")
}

writeChangeoDb(results, file=output)

