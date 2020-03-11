# paper: http://dx.doi.org/10.1093/nar/gkv1198
# docs from supplementary info:

# issues
#  - can't pass in your own tree -- have to let it run either phylip dnapars or phylip neighbor. Related things:
#    - it assumes the internal nodes (that phylip adds) have phylip-convention names, i.e. str(some_integer), and if they're not, it breaks
#    - if you pass it (i.e. if you observe) internal nodes, phylip shunts them to zero-lengh branches hanging off of newly-inferred internal nodes. So you then have to try to figure out which nodes/leaves you can collapse
#  - it counts offspring all the way down the tree, i.e. ancestors get credit for fitness improvements in their offspring
#    - this is "handled" by the "affected by descendents" flag, but that doesn't solve the problem, it just alerts you to lineages where there's more likely to be a problem
#    - it would be much better to use an lbi-style weighting that decreases contributions as they get further away
#  - ignores pairs of sibling edges in which both edges have mutations, i.e. it can throw out a large fraction of mutations in cases where most branches have mutations
#  - ignores branches with no siblings
#  - ignores branches with multiple mutations
# notes
#  - number of offspring for a node is set as the number of edges in the entire subtree below that node

## Computing LONR scores in lineage trees

## General
##   The program is written in R.
##   This analysis is divided into two parts:
##     First, a lineage tree is built with sequences provided in and aligned FASTA format.
##       There is an option to provide an outgroup sequence in order to detect more accurately the root sequence.
##       For less than 100 sequences, the Maximum parsimony method is used (dnapars Phylip program).
##       Otherwise, the Neighbor-joining method is used (neighbor  Phylip program).
##       The Neighbor-joining method uses a distance matrix as input, which is computed using the dnadist Phylip program.
##       The internal sequences are reconstructed using the Fitch algorithm.
##       The tree is then divided into subtrees in order to ignore mutations occurring in very distant sequences.
##       Thus the user can specify the cut-off (in nucleotides) to cut long branches.
##       The default is 10 mutations.
##       The following LONR analysis is then performed for each subtree separately.
##     Second, mutations are detected between each pair of father and son sequences (at the nucleotide level).
##       The LONR score is calculated, for each mutation, as the log of the ratio between the sub-tree size of the son in which the mutation occurred and the sub-tree size of the son in which no mutation occurred at this position.

## How to run
##   This program uses the Phylip-3.695 package. It is expecting it to be in the same folder as this script. The path to the program can be modified in the beginning of the script.
##   The package Biostrings and seqinr need to be installed.

## Main function  - compute.LONR()

## Input
##   1) in.dir – aligned FASTA file input directory
##   2) out.dir – output directory for tree files and LONR results
##   3) file – FASTA file name
##   4) outgroup (optional) – outgroup sequence name. If
##   5) cutoff  - branches with more mutations than specified nu this cutoff will be trimmed, resulting in several subtrees (default  - 10 )

## Output
##  The following directories are created in out.dir:
##   Tree directory –
##     1) filename.fasta -modified FASTA file if gaps were removed.
##     2) filename.dis – (only for neighbor joining trees) contains distance matrix created by the dnadist program in the Phylip package.
##     3) filename.phy – original sequences in alignment format.
##     4) filename_edges.tab – tab-delimited file containing the tree edges, their weights and the distance in nucleotides between each two nodes.
##     5) filename_names.tab – tab-delimited file matching between original sequence names and temporary names.
##     6) filename_out.txt – Phylip output tree file
##     7) filename_tree.txt – Phylip output tree file in Newick format

## LONR directory –
##  1) filename_lonr.csv – comma-separated file containing lonr results as followed:
##    i.   mutation – mutated nucleotides (e.g. AC means A C)
##    ii.  LONR – log(size of mutated sub-tree/size of un-mutated sub-tree)
##    iii. mutation.type – (S) Silent or (Replacement)
##    iv.  position – in nucleotides, according to output FASTA file (see above)
##    v.   father – sequence name from which occurred mutation
##    vi.  son – sequence name to which occurred mutation
##    vii. flag – True if LONR score an internal node is affected by mutations occurring in its descendants
## Clarifications
##  1) The input sequences must already be aligned. If there are gaps, the consensus sequence of all the sequences is computed, and positions containing gaps are removed from all the sequences. The output FASTA file is created after this step.
##  2) The dnapars program may create trees which are not completely binary. Thus, internal nodes which have more than two children are fixed by created an identical new child, which will receive the extra children.
##     This also happens is case no outgroup is provided and the root has three children (both in dnapars and neighbor programs).
##  3) If the nucleotide sequence lengths are not a multiple of three and mutations occurred in the last nucleotides, these mutations are ignored since they cannot be typed (not a full codon).

# imports
suppressPackageStartupMessages(require(seqinr, quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(require(Biostrings, quietly=TRUE, warn.conflicts=FALSE))

MIN.SEQ <- 3
MAX.SEQ <- 7000

# ----------------------------------------------------------------------------------------
# empty ones are all set by the calling python script after importing this code
G.phy.infname = 'inseqs.phy'
G.dnadist.fname = 'dnadist.dis'
G.phy.outfname = ''
G.phy.treefname = ''
G.outseqs.fname = ''
G.edgefname = ''
G.names.fname = ''
G.lonrfname = ''

# Create fasta file from a data.frame of sequences and headers
#          out.dir - output directory
#          nameSeq.df  - data.frame of headers and sequences
write.FASTA <- function(out.dir, nameSeq.df) {
  sequences <- nameSeq.df$seq
  names(sequences) <- nameSeq.df$head
  writeXStringSet(DNAStringSet(sequences), file=paste0(out.dir, G.outseqs.fname), width=1000)
}

# Change sequence names
#          fasta.df - data frame of sequences and headers
#          outgroup - outgroup sequence name
# Returns:  data frame of sequences and headers with extra column of new headers
change.names <- function(fasta.df, outgroup = NULL){

  if (!is.null(outgroup)){
    outgroup.ind <- which(fasta.df$head == outgroup)
    n.seq <- nrow(fasta.df)
    # change sequence names to 'L_#'
    ind <- 1:n.seq
    ind <- ind[-outgroup.ind]

    # original method, fails if outgroup is the last sequence in the file (it complains that the first line only adds n.seq-1 entries, whereas the rest of the data frame has n.seq entries [or something like that]):
    ## fasta.df$head2[ind] <- sapply(1:(n.seq-1), function(x) { paste0('L', x)})
    ## fasta.df$head2[outgroup.ind] <- fasta.df$head[outgroup.ind]

    # this is weird and complicated so that it duplicates the original method ^
    fasta.df$head2[1:n.seq] <- sapply(1:n.seq, function(x) { if(x == outgroup.ind) return(fasta.df$head[outgroup.ind]); if(x > outgroup.ind) x = x-1; return(paste0('L', x))})
  }else{
    n.seq <- nrow(fasta.df)
    # change sequence names to 'L_#'
    ind <- 1:n.seq
    fasta.df$head2[ind] <- sapply(1:n.seq, function(x) { paste0('L', x)})
  }
  return(fasta.df)
}

# Create phylip input alignment file
# Arguments:   fasta.df - data frame of sequences and headers
fasta2phylip <- function(fasta.df, workdir){

  phy.df <- rbind(data.frame(head=sprintf('%-9s', nrow(fasta.df)), seq=nchar(fasta.df$seq[1]), stringsAsFactors=F),
                  data.frame(head=sprintf('%-9s', fasta.df$head2), seq=fasta.df$seq, stringsAsFactors=F))

  write.table(phy.df, file=paste0(workdir, G.phy.infname), quote=F, sep=' ', col.names=F, row.names=F)
}

# Run dnapars (maximum parsimony)
#           outgroup.ind - outgroup index in file (optional)
run.dnapars <- function(workdir, outgroup.ind){
  print('  building trees with maximum parsimony')

  curr.dir <- getwd()
  setwd(workdir)

  # options from dnapars program
  # index of outgroup - last in data frame
  if (length(outgroup.ind) != 0)
    pars.options <- c(G.phy.infname, 'O', outgroup.ind, 'V', '1', '5', '.', '2', 'Y')
  else
    pars.options <- c(G.phy.infname, 'V', '1', '5', '.', '2', 'Y')

  # run dnapars
  system2('phylip', args='dnapars', input=pars.options, stdout=NULL)  # not sure why this was just calling 'dnapars'? my phylip install doesn't seem to have put dnapars in my path, but 'phylip dnapars' seems to work ok

  # move .phy file and output tree files
  file.rename(from=paste0(workdir, 'outfile'), to=paste0(workdir, G.phy.outfname))
  file.rename(from=paste0(workdir, 'outtree'), to=paste0(workdir, G.phy.treefname))

  setwd(curr.dir)
}

# Order sequence from root to leaves
#
# Arguments:   root - root sequence name
#           nameSeq.df - data frame with headers and sequences
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
#           outgroup  - outgroup sequence name (optional)
#
# Returns:  nameSeq.df - data frame with headers and sequences, ordered by distance from root
order.nameSeq <- function(root, nameSeq.df, edge.df, outgroup){

  n.seq <- nrow(nameSeq.df)
  nameSeq.df2 <- nameSeq.df[which(nameSeq.df[,'head']==root),] # root
  for(i in 1:n.seq){ #1 - root
    sons <- edge.df[which(edge.df[,'from']==nameSeq.df2[i,'head']),2]
    if(length(sons)>0){
      outgroup.ind <- which(sons==outgroup)
      if(length(outgroup.ind)!=0) { sons <- sons[-outgroup.ind] }
      for(j in 1:length(sons))
        nameSeq.df2 <- rbind(nameSeq.df2, nameSeq.df[which(nameSeq.df[,'head']==sons[j]),])
    }
  }
  if (!is.null(outgroup)) # last - outgroup
    nameSeq.df2 <- rbind(nameSeq.df[which(nameSeq.df[,'head']==outgroup),], nameSeq.df2) # add outgroup first
  row.names(nameSeq.df2)<-NULL

  return(nameSeq.df2)
}


# Parse dnapars output tree file, get internal sequences and tree structure
#
# Returns:  list of:
#           nameSeq.df - data frame with headers and sequences, ordered by distance from root
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
parse.dnapars <- function(workdir, outgroup = NULL){

  # read output tree file
  out.tree <- scan(paste0(workdir, G.phy.outfname), what='character',sep='\n',
                   blank.lines.skip=F, strip.white=F)
  # check if tree was build
  if (any(grepl('-1 trees in all found', out.tree))) { return(NULL) }

  # get internal sequences
  seq.start <- min(grep('From\\s+To\\s+Any Steps\\?\\s+State at upper node', out.tree, perl=T, fixed=F))
  seq.empty <- grep('^\\s*$', out.tree[seq.start:length(out.tree)], perl=T, fixed=F)
  seq.len <- seq.empty[min(which(seq.empty[-1] == (seq.empty[-length(seq.empty)] + 1)))]
  seq.block <- paste(out.tree[(seq.start + 2):(seq.start + seq.len - 2)], collapse='\n')
  seq.df <- read.table(textConnection(seq.block), as.is=T, fill=T, blank.lines.skip=T,  row.names = NULL, header=F, stringsAsFactors=F)

  # fix root lines and remove empty rows
  fix.row <- which(seq.df[,3]!="yes" & seq.df[,3]!="no" & seq.df[,3]!="maybe")
  if (!is.null(outgroup))
    seq.df[fix.row, ] <- cbind(seq.df[fix.row, 1], seq.df[fix.row, 2],'no', seq.df[fix.row, 3:6], stringsAsFactors=F)
  else
    seq.df[fix.row, ] <- cbind(seq.df[fix.row, 1], seq.df[fix.row, 1],'no', seq.df[fix.row, 2:5], stringsAsFactors=F)

  # save full sequences as a data frame
  names <- unique(seq.df[, 2])
  seq <- sapply(names, function(x) { paste(t(as.matrix(seq.df[seq.df[, 2] == x, -c(1:3)])), collapse='') })
  nameSeq.df <- data.frame(head=names, seq=seq, stringsAsFactors=F, row.names=NULL)

  # get tree structure
  edge.start <- min(grep('between\\s+and\\s+length', out.tree, perl=T, fixed=F))
  edge.len <- min(grep('^\\s*$', out.tree[edge.start:length(out.tree)], perl=T, fixed=F))
  edge.block <- paste(out.tree[(edge.start + 2):(edge.start + edge.len - 2)], collapse='\n')
  edge.df <- read.table(textConnection(edge.block), col.names=c('from', 'to', 'weight'), as.is=T, stringsAsFactors=F)

  # order sequences by distance from root
  root <- unique(edge.df$from)[!(unique(edge.df$from) %in% edge.df$to)]

  #root <- seq.df[which(seq.df[,1]=='root')[1],2]
  nameSeq.df <- order.nameSeq(root, nameSeq.df, edge.df, outgroup)

  return(list(nameSeq.df, edge.df))
}

# Fix ambiguous nucleotides in internal sequences, according to IUPAC code
#
# Arguments:   nameSeq.df - data frame with headers and sequences, order by distance from root
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
#           outgroup - outgroup sequence name (optional)
#
# Returns:  - nameSeq.df - data frame with fixed sequences
fix.internal <- function(nameSeq.df, edge.df, outgroup = NULL){

  nucleotides <- c('A', 'C','G','T','-')

  # from leaves to root (except outgroup)
  # fix only internal nodes!
  if (!is.null(outgroup))
    intern.seq <- rev(setdiff(2:nrow(nameSeq.df),grep('^L', nameSeq.df[,'head'], perl=T, fixed=F)))
  else
    intern.seq <- rev(setdiff(1:nrow(nameSeq.df),grep('^L', nameSeq.df[,'head'], perl=T, fixed=F)))

  for(i in intern.seq){
    # get ambiguous nucleotides
    curr.seq <- nameSeq.df[i,'seq']
    amb.pos <- setdiff(1:nchar(curr.seq),unlist(lapply(nucleotides, function(x) {gregexpr(pattern =x,curr.seq)})))
    if(length(amb.pos) > 0){
      # get son sequences
      sons.seq <- sapply(edge.df[which(edge.df[,'from']==nameSeq.df[i,'head']),'to'], function(x) {nameSeq.df[which(nameSeq.df[,'head']==x),'seq']})
      # do not fix outgroup sequence
      if (!is.null(outgroup)){
        outgroup.ind <- which(names(sons.seq)==outgroup)
        if(length(outgroup.ind)!=0)
          sons.seq <- sons.seq[-outgroup.ind]
      }

      for(j in amb.pos){
        amb.nuc <- substr(curr.seq,j,j)
        sons.nuc <- substr(sons.seq,j,j)
        # if deletion
        if(amb.nuc == 'O')
          curr.seq <- paste0(substr(curr.seq,1, j-1), '-', substr(curr.seq,j+1, nchar(curr.seq)))
        else{
          # check if at least one of the sons has one match for ambiguous letter
          if(is.element(amb.nuc, c('?', 'X')))
            matches <- which(is.element(sons.nuc, nucleotides))
          else
            matches <- which(is.element(sons.nuc, toupper(amb(amb.nuc, forceToLower = TRUE))))
          if(length(matches)==1) # if only one son had a match
            curr.seq <- paste0(substr(curr.seq,1, j-1), sons.nuc[matches], substr(curr.seq,j+1, nchar(curr.seq)))
          else if(length(matches)==0){ # choose randomly from son nucleotides
            rand <- sample(1:length(sons.nuc), 1, replace=T)
            curr.seq <- paste0(substr(curr.seq,1, j-1), sons.nuc[rand], substr(curr.seq,j+1, nchar(curr.seq)))
          }else{ # choose randomly from son nucleotides that matched
            rand <- sample(1:length(matches), 1, replace=T)
            curr.seq <- paste0(substr(curr.seq,1, j-1), sons.nuc[matches[rand]], substr(curr.seq,j+1, nchar(curr.seq)))
          }
        }
      }
      nameSeq.df[i,'seq'] <- curr.seq
    }
  }

  return(nameSeq.df)
}

# Compute distance matrix (Kimura model) for neighbor joining using dnadist (Phylip package)
#
# Arguments:
#          workdir - phylip directory
run.dnadist <- function(workdir){

  # options from dnadist program
  dnadist.options <- c(G.phy.infname, 'D', '2', 'Y')

  # run dnadist
  system2('phylip', args='dnadist', input=dnadist.options,, stdout=NULL)

  # move .phy file and output tree files
  file.rename(from=paste0(workdir, 'outfile'), to=paste0(workdir, G.dnadist.fname))
}


# Run neighbor (Neighbor joining)
#          outgroup.ind - outgroup index  in file (optional)
run.neighbor <- function(workdir, outgroup.ind){
  print('  building trees with neighbor joining')
  curr.dir <- getwd()
  setwd(workdir)

  # run dnadist to create distance matrix
  run.dnadist(workdir)

  # options for neighbor program
  # index of outgroup - last data frame
  if (length(outgroup.ind) != 0 )
    neigh.options <- c(G.dnadist.fname, 'O', outgroup.ind, '2', 'Y')
  else
    neigh.options <- c(G.dnadist.fname, '2', 'Y')

  # run neighbor
  system2('phylip', args='neighbor', input=neigh.options, stdout=NULL)

  # move .phy, .dis and output tree files
  file.rename(from=paste0(workdir, 'outfile'), to=paste0(workdir, G.phy.outfname))
  file.rename(from=paste0(workdir, 'outtree'), to=paste0(workdir, G.phy.treefname))
  file.remove(paste0(workdir, G.dnadist.fname))

  setwd(curr.dir)
}

# Parse neighbor output tree file, get internal sequences and tree structure
#           fasta.df - data frame of input sequences and headers
#           outgroup - outgroup sequence name (optional)
# Returns:  - list of:
#           nameSeq.df - data frame with headers and sequences
#           edge.df - data frame with columns - parent(from), child(to), distance(weight)
parse.neighbor <- function(workdir, fasta.df, outgroup = NULL){

  # read output tree file
  out.tree <- scan(paste0(workdir, G.phy.outfname), what='character',sep='\n',
                   blank.lines.skip=F, strip.white=F)

  # check if tree was build
  if (any(grepl('-1 trees in all found', out.tree))) { return(NULL) }

  # get tree structure
  edge.start <- min(grep('Between\\s+And\\s+Length', out.tree, perl=T, fixed=F))
  edge.len <- min(grep('^\\s*$', out.tree[edge.start:length(out.tree)], perl=T, fixed=F))
  edge.block <- paste(out.tree[(edge.start + 2):(edge.start + edge.len - 2)], collapse='\n')
  edge.df <- read.table(textConnection(edge.block), col.names=c('from', 'to', 'weight'), as.is=T)

  # create nameSeq.df from input sequence only (for now)
  nameSeq.df <- data.frame(head = fasta.df[,'head2'], seq = fasta.df[,'seq'], stringsAsFactors=F)

  return(list(nameSeq.df, edge.df))
}

# Traverse tree in postorder (step 2 in Fitch's algorithm)
#
# Arguments:  father/ sons - father and its 2 son names
#          edge.df - data frame with columns - parent(from), child(to), distance(weight)
#          nameSeq.list - list of sequences with headers
#
# Returns: modified nameSeq.list - list of sequences with headers
traverse.up <- function(father, sons, edge.df, nameSeq.list){

  sons.seq <- sapply(sons,function(x) NULL)

  # check if each son sequence a leaf or is already reconstructed
  for(i in sons){
    if(!is.element(i, names(nameSeq.list)))
      nameSeq.list <-  traverse.up(i, edge.df[edge.df[,'from']==i,'to'], edge.df, nameSeq.list)
    sons.seq[i] <- nameSeq.list[i]
  }
  seq.len <- length(sons.seq[[1]])

  father.seq <- character(seq.len)
  for(i in 1:seq.len){ # for each position in sequence
    curr.nuc <- sapply(sons.seq, function(x)  x[[i]])
    nuc <- Reduce(intersect, strsplit(curr.nuc,""))
    if(length(nuc)>0) # save intersection
      father.seq[i] <- paste(nuc, collapse='')
    else # save union
      father.seq[i] <- paste( Reduce(union, strsplit(curr.nuc,'')),collapse="")
  }
  nameSeq.list[[father]]<-father.seq # save modified sequence

  return(nameSeq.list)
}

# Recursively traverse tree in preorder (step 1 in Fitch's algorithm)
#
# Arguments:  father/ son - father and son name
#          edge.df - data frame with columns - parent(from), child(to), distance(weight)
#          nameSeq.list - list of sequences with headers
#          outgroup - outgroup sequence name (optional)
#
# Returns: modified nameSeq.list - list of sequences with headers
traverse.down <- function(father, son, edge.df, nameSeq.list, outgroup = NULL){

  son.seq <- nameSeq.list[[son]]
  seq.len <- length(son.seq)

  if(is.null(father)){ # root
    for(i in 1:seq.len){
      if(nchar(son.seq[i])>1){ # more than 1 option in intersection - choose randomly
        #rand <- floor(runif(1, min=1, max=nchar(son.seq[i])+1))
        rand <- sample(1:nchar(son.seq[i]), 1, replace=T)
        son.seq[i] <- substr(son.seq[i],rand, rand)
      }
    }
  }else{
    father.seq <- nameSeq.list[[father]]
    for(i in 1:seq.len){
      if(nchar(son.seq[i])>1){
        nuc <- Reduce(intersect, strsplit(c(father.seq[i],son.seq[i]),''))
        # if only one nucleotide in intersection - keep it
        if(length(nuc)==0){  # no intersection - choose randomly from son's options
          rand <- sample(1:nchar(son.seq[i]), 1, replace=T)
          son.seq[i] <- substr(son.seq[i],rand, rand)
        }else # if intersection is not empty, but son has more than 1 option
          son.seq[i] <- nuc
      }
    }
  }
  nameSeq.list[[son]] <- son.seq # save modified sequence

  # recursive call for each son's sons
  sons <- edge.df[edge.df[,'from']==son,'to']
  if(length(sons) > 0){
    # remove outgroup from root's children
    if (!is.null(outgroup)){
      outgroup.ind <- which(sons == outgroup)
      if(length(outgroup.ind)!=0) { sons <- sons[-outgroup.ind] }
    }
    nameSeq.list <- traverse.down(son, sons[1], edge.df, nameSeq.list, outgroup)
    nameSeq.list <- traverse.down(son, sons[2], edge.df, nameSeq.list, outgroup)
  }
  return(nameSeq.list)
}

# Reconstruct internal sequences with Fitch algorithm (for neighbor trees)
#
# Arguments:  nameSeq.df - data frame with headers and input sequences only
#          edge.df - data frame with columns - parent(from), child(to), distance(weight)
#          outgroup - outgroup sequence name (optional)
#
# Returns: nameSeq.df - data frame with headers and all sequences, ordered by distance from root
get.internal <- function(nameSeq.df, edge.df, outgroup = NULL){

  # find root
  root <- as.character(unique(edge.df[!(edge.df[, 1] %in% edge.df[, 'to']), 1]))
  root.sons <- edge.df[edge.df[,'from']==root,'to']

  # remove outgroup from root's children
  if (!is.null(outgroup)){
    outgroup.ind <- which(root.sons == outgroup)
    if(length(outgroup.ind)!=0) { root.sons <- root.sons[-outgroup.ind] }
  }

  # convert sequences to lists
  nameSeq.list <- sapply(nameSeq.df[,'seq'], function(x) {strsplit(x, "", fixed=FALSE)})
  names(nameSeq.list) <-  nameSeq.df[,'head']

  # Step 1 - preorder on tree - get intersection, otherwise - get union
  nameSeq.list <- traverse.up(root, root.sons, edge.df, nameSeq.list)

  # Step 2 - postorder on tree - get intersection, otherwise choose randomly
  nameSeq.list <- traverse.down(NULL, root, edge.df, nameSeq.list, outgroup)

  # convert sequence list to dataframe
  nameSeq.list2 <-lapply(nameSeq.list,function(x) {paste(x, collapse="")})
  nameSeq.df <-data.frame(head = names(nameSeq.list2), stringsAsFactors=F)
  nameSeq.df$seq <- unlist(nameSeq.list2)

  # order sequences in nameSeq.df by distance from root
  nameSeq.df <- order.nameSeq(root, nameSeq.df, edge.df, outgroup)

  return(nameSeq.df)
}

# Match old input sequence names with new in final output sequence data frame
#
# Arguments:  nameSeq.df - data frame with headers and all sequences
#          fasta.df - data frame with new and old headers and input sequences only
#
# Returns: nameSeq.df - data frame with headers and sequences, with old and new names
match.names <- function(fasta.df, nameSeq.df){

  nameSeq.df$head2 <- rep("-", nrow(nameSeq.df))
  for(i in 1:nrow(fasta.df))
    nameSeq.df[which(fasta.df[i,'head2']==nameSeq.df[,'head']),'head2'] <- fasta.df[i,'head']

  return(nameSeq.df)
}

# Compute edge lengths (distance in nucleotides)
#
# Arguments:  nameSeq.df - data frame with headers and sequences
#          edge.df - data frame with columns - parent(from), child(to), edge weight (weight)
#
# Returns: edge.df - data frame with columns - parent(from), child(to), edge weight (weight) , edge length (distance in nt)
compute.edge <- function(nameSeq.df, edge.df){

  n.edge <- nrow(edge.df)
  edge.df$distance <- rep(0, n.edge)
  for(i in 1:n.edge){
    from.seq <- unlist(strsplit(nameSeq.df[nameSeq.df[,'head']==as.character(edge.df[i,'from']),'seq'],''))
    to.seq <- unlist(strsplit(nameSeq.df[nameSeq.df[,'head']==as.character(edge.df[i,'to']),'seq'],''))
    edge.df[i,'distance'] <- length(which(from.seq!=to.seq))
  }

  return(edge.df)
}

# Convert parsimony tree into binary tree
#
# Arguments:  nameSeq.df - data frame with headers and sequences
#          edge.df - data frame with columns - parent(from), child(to), edge weight (weight)
#          outgroup name (default - NULL)
#
# Returns: modified edge.df and nameSeq.df
convert.to.binary <- function(nameSeq.df, edge.df, outgroup = NULL){

  # find nodes with more than 2 sons
  son.count <- table(edge.df$from)
  nodes <- names(son.count)[son.count == 3]

  # get name of last internal node created by program
  ind <- max(unique(edge.df$from)) + 1

  for (i in nodes){
    sons <- edge.df[edge.df$from==i, 'to']
    if (!is.null(outgroup)){
      if (outgroup %in% sons) # if outgroup is one of the root's sons - ok
        next
    }

    # create new node named with current ind and attach 2 sons
    edge.df[edge.df$to==sons[1], 'from'] <- ind
    edge.df[edge.df$to==sons[2], 'from'] <- ind
    # add new node has son of current node
    edge.df <- rbind(edge.df, data.frame(from = i, to = ind, weight=0))
    # add new node sequence - same as current node
    nameSeq.df <- rbind(nameSeq.df, data.frame(head=ind, seq=nameSeq.df[nameSeq.df$head==i, 'seq']))
    ind <- ind + 1
  }

  return(list(nameSeq.df, edge.df))

}

# Build phylogenetic tree with maximum parsimony or neighbor joining
# Arguments:
#          fasta.df - data frame with input sequences and headers
#          outgroup - outgroup sequence name (optional)
# Returns: nameSeq.list - list of sequences with headers
#          edge.df - data frame with columns - parent(from), child(to), weight, distance(nt)
build.trees <- function(method, fasta.df, workdir, outgroup=NULL, existing.edgefile=NULL, existing.node.seqfile=NULL){

  dir.create(workdir, recursive = T, showWarnings = FALSE)

  # change input sequence name to shorter ones
  fasta.df <- change.names(fasta.df, outgroup) # CHECK CP NUM!!!!!!!

  # convert sequences to alignment format for phylip programs
  fasta2phylip(fasta.df, workdir)

  # find outgroup index in file
  outgroup.ind <- which(fasta.df$head == outgroup)

  if(!is.null(existing.edgefile)) {
    stop('doesn\'t work yet (and probably a waste of time to fix it -- reimplement this stupid shit in python)')
    edge.df <- read.csv(file=existing.edgefile, header=TRUE, sep=',', colClasses=c('character', 'character', 'numeric'))
    nameSeq.df <- read.csv(file=existing.node.seqfile, header=TRUE, sep=',', colClasses=c('character', 'character'))
    for(irow in 1:length(nameSeq.df$head)) {
        input.name <- nameSeq.df$head[[irow]]
        for(jrow in 1:length(fasta.df$head)) {
            if(fasta.df$head[[jrow]] == input.name) {  # see if this sequence is in the df with name translations (i.e. it's an input/leaf sequence)
                new.name <- fasta.df$head2[[jrow]]
                nameSeq.df$head[[irow]] <- new.name  # and if it is, replace the input  name with the new name
                if(toString(nameSeq.df$seq[[irow]]) != toupper(toString(fasta.df$seq[[jrow]])))  # make sure the sequences are the same (I don't really understand the need for toString(), but it prints some extra "Level" stuff if I don't have it, so...)
                    stop(paste0('sequences don\'t match for ', input.name, ':\n', toString(nameSeq.df$seq[[irow]]), '\n', toupper(toString(fasta.df$seq[[jrow]]))))
                for(krow in 1:length(edge.df$from)) {  # also have to change the names in the edge df
                    if(edge.df$from[[krow]] == input.name)  # this doesn't actually happen, since they only rename the leaves
                        edge.df$from[[krow]] <- new.name
                    if(edge.df$to[[krow]] == input.name)
                        edge.df$to[[krow]] <- new.name
                }
                break
            }
        }
    }
    ## # we need to convert the internal node names to integers here, since convert.to.binary() tries sorts them, which fails if they're strings. But this doesn't freaking work, because the class of the column is set the 'character' (above). fuck you, R
    ## print(nameSeq.df$head)
    ## ftmp <- function(x) { if(substr(x, 1, 1) == 'L') return(x); if(x == outgroup) return(x); return(strtoi(x)); }
    ## edge.df$from <- sapply(edge.df$from, ftmp)
    ## edge.df$to <- sapply(edge.df$to, ftmp)
    ## nameSeq.df$head <- sapply(nameSeq.df$head, ftmp)
    ## print(nameSeq.df$head)
    ## stop('x')
  } else if(method == 'dnapars') {  # maximum parsimony
    run.dnapars(workdir, outgroup.ind)

    # parse dnapars output
    tmp <- parse.dnapars(workdir, outgroup)
    if(is.null(tmp)) { return(F) } # if tree building failed
    nameSeq.df <- tmp[[1]]
    edge.df <- tmp[[2]]
    nameSeq.df <- fix.internal(nameSeq.df, edge.df, outgroup)
  } else if(method == 'neighbor') {  # neighbor joining
    run.neighbor(workdir, outgroup.ind)

    # parse neighbor output
    tmp <- parse.neighbor(workdir, fasta.df, outgroup)
    if(is.null(tmp)) { return(F) } # if tree building failed
    nameSeq.df <- tmp[[1]]
    edge.df <- tmp[[2]]
    nameSeq.df <- get.internal(nameSeq.df, edge.df, outgroup)
  } else {
    stop(paste0('unhandled method: ', method))
  }

  # convert tree to binary
  tmp <- convert.to.binary(nameSeq.df, edge.df, outgroup)
  nameSeq.df <- tmp[[1]]
  edge.df <- tmp[[2]]

  # compute edge lengths
  edge.df <- compute.edge(nameSeq.df, edge.df)

  # retrieve old names for input sequences
  nameSeq.df <- match.names(fasta.df, nameSeq.df)

  # save output files
  # save sequence as FASTA file
  write.FASTA(workdir, nameSeq.df)
  # save Fome/To/distances table (edges)
  write.table(edge.df, file=paste0(workdir, G.edgefname), quote=F, sep='\t', col.names=T, row.names=F)
  # save old and new sequence names
  write.table(nameSeq.df[,c('head','head2')], file=paste0(workdir, G.names.fname), quote=F, sep='\t', col.names=T, row.names=F)

  return(list(nameSeq.df, edge.df))
}

# Calculate subtree sizes for each internal node (recursive function)
#
# Arguments:  node - current node for which the number of offsprings
#          n.offspring - data frame with size of keft and right subtrees for each sequence
#          edge.df - data frame with columns - parent(from), child(to), weight, distance(nt)
#
# Returns: updated n.offspring
get.number.offsprings <- function(node, n.offspring, edge.df){

  # post order on tree
  sons <- edge.df[edge.df[, 'from'] == node, 2]

  # if node is a leaf - set number of offsprings to zero
  if (length(sons) == 0)
    return (n.offspring)

  # sort sons alphabetically - so that the first is always the left son
  sons <- sort(sons)

  node.ind <- which(n.offspring$node == node)

  # left subtree
  n.offspring <- get.number.offsprings(sons[1], n.offspring, edge.df)
  n.offspring[node.ind, 'left'] <- n.offspring[n.offspring$node == sons[1], 'left'] +
    n.offspring[n.offspring$node == sons[1], 'right'] + 1

  if (length(sons) == 1 ) # the other son was cut
    n.offspring[node.ind, 'right'] <- 0
  else{
    # right subtree
    n.offspring <- get.number.offsprings(sons[2], n.offspring, edge.df)
    n.offspring[node.ind, 'right'] <- n.offspring[n.offspring$node == sons[2], 'left'] +
      n.offspring[n.offspring$node == sons[2], 'right'] + 1
  }

  return (n.offspring)
}

# Analyze mutations between a pair of father-son and compute LONR scores for each mutation
#
# Arguments:  n.offspring - data frame with size of left and right subtrees for each sequence
#          father - father name
#          father.char - father nucleotide sequence split into characters
#          father.aa.char - father amino acid sequence split into characters
#          son - son name
#          son.side - left or right son of father (alphabetically)
#          son.char - son nucleotide sequence split into characters
#          son.aa.char - son amino acid sequence split into characters
#          mut.pos - mutation positions (nucleotides)
#          mutations - mutation table
#          mutations.ind - index of last mutation inserted in mutations table
#
# Returns: list of
#           - updated mutations and mutations.ind
analyze.mutations <- function(n.offspring, father, father.char, father.aa.char, son, son.side, son.char, son.aa.char, mut.pos, mutations, mutations.ind){

  flag <- F
  father.ind <- which(n.offspring$node==father)

  # if at least one son was trimmed - do not analyze mutations
  if (n.offspring[father.ind,'left']==0 | n.offspring[father.ind,'right']==0)
    return(list(mutations, mutations.ind))

  # compute LONR score
  # log of subtree size where mutation occurred divided by subtree size where no mutation occurred
  if (son.side == 'left')
    lonr <- log(n.offspring[father.ind,'left']/n.offspring[father.ind,'right'])
  else
    lonr <- log(n.offspring[father.ind,'right']/n.offspring[father.ind,'left'])

  # get son's subtree sizes
  son.left <- n.offspring[n.offspring$node==son, 'left']
  son.right <- n.offspring[n.offspring$node==son, 'right']
  if (son.left != 0 & son.right != 0){
    if (son.left > 2*son.right | 2*son.left < son.right){
      if (son.side == 'left')
        new.lonr <- log((2*min(son.left, son.right))/n.offspring[father.ind,'right'])
      else
        new.lonr <- log((2*min(son.left, son.right))/n.offspring[father.ind,'left'])
      if (lonr/new.lonr < 0) # sign changed - add flag
        flag = T
    }
  }

  for (j in mut.pos){
    # mutation (nt)
    mutations[mutations.ind,'mutation'] <- paste0(father.char[j],son.char[j])
    # LONR
    mutations[mutations.ind,'LONR'] <- lonr
    # mutation type
    # if mutation occurred in last nucleotides which are not a full codon - do not analyse
    if( length(father.aa.char) >= ceiling(j/3) ){
      if (father.aa.char[ceiling(j/3)] != son.aa.char[ceiling(j/3)])
        mutations[mutations.ind,'mutation.type'] <- 'R'
      else
        mutations[mutations.ind,'mutation.type'] <- 'S'
    }else
      next
    # position (nt)
    mutations[mutations.ind,'position'] <- j

    # save father and son names
    mutations[mutations.ind,'father'] <- father
    mutations[mutations.ind,'son'] <- son
    mutations[mutations.ind,'flag'] <- flag

    mutations.ind <- mutations.ind + 1

    # increase table size if needed
    if (mutations.ind > nrow(mutations))
      mutations <- rbind(mutations, data.frame(mutation = rep('', 1000), LONR = 0, mutation.type = '', position=0, father='', son='', flag=F, stringsAsFactors=F))
  }

  return(list(mutations, mutations.ind))
}


# Get mutations between each father-son pair (recursively)
#
# Arguments:  node - current node for which the number of offsprings
#          edge.df - data frame with columns - parent(from), child(to), weight, distance(nt)
#          n.offspring - data frame with size of left and right subtrees for each sequence
#          mutations - mutation table
#          mutations.ind - index of last muation inserted in mutations table
#
# Returns: updated n.offspring
get.mutations <- function(node, edge.df, nameSeq.df, n.offspring, mutations, mutations.ind){

  sons <- edge.df[edge.df[, 1] == node, 2]

  # if node is a leaf -
  if (length(sons) == 0)
    return (list(mutations, mutations.ind))

  sons <- sort(sons)

  # get father sequence, translate and split into characters
  father.seq <- nameSeq.df[nameSeq.df$head==node, 'seq']
  father.char <- s2c(father.seq)
  father.aa.char <- seqinr::translate(unlist(strsplit(tolower(substr(father.seq, 1, nchar(father.seq)-(nchar(father.seq)%%3))), "")),
                                      numcode = 1, NAstring = "X", ambiguous = FALSE)

  # get son sequences, translate and split into characters
  son1.seq <- nameSeq.df[nameSeq.df$head==sons[1], 'seq']
  son1.char <- s2c(son1.seq)
  son1.aa.char <- seqinr::translate(unlist(strsplit(tolower(substr(son1.seq, 1, nchar(son1.seq)-(nchar(son1.seq)%%3))), "")),
                                    numcode = 1, NAstring = "X", ambiguous = FALSE)
  # get mutation positions
  mut.pos.son1 <- which(father.char!=son1.char)

  if (length(sons) > 1){
    # get son sequences, translate and split into characters
    son2.seq <- nameSeq.df[nameSeq.df$head==sons[2], 'seq']
    son2.char <- s2c(son2.seq)
    son2.aa.char <- seqinr::translate(unlist(strsplit(tolower(substr(son2.seq, 1, nchar(son2.seq)-(nchar(son2.seq)%%3))), "")),
                                      numcode = 1, NAstring = "X", ambiguous = FALSE)
    # get mutation positions
    mut.pos.son2 <- which(father.char!=son2.char)
    # do not analyze mutation in both sons
    mut.pos.son1 <- setdiff(mut.pos.son1,intersect(mut.pos.son1, mut.pos.son2))
    mut.pos.son2 <- setdiff(mut.pos.son2,intersect(mut.pos.son1, mut.pos.son2))
  }
  # analyze mutation in first son
  res<- analyze.mutations(n.offspring, node, father.char, father.aa.char, sons[1], 'left', son1.char, son1.aa.char, mut.pos.son1, mutations, mutations.ind)
  mutations <- res[[1]]
  mutations.ind <- res[[2]]
  # recursive call with sons
  res <- get.mutations(sons[1], edge.df, nameSeq.df, n.offspring, mutations, mutations.ind)
  mutations <- res[[1]]
  mutations.ind <- res[[2]]

  if (length(sons) > 1){
    # analyze mutation in second son
    res<- analyze.mutations(n.offspring, node, father.char, father.aa.char, sons[2], 'right', son2.char, son2.aa.char, mut.pos.son2, mutations, mutations.ind)
    mutations <- res[[1]]
    mutations.ind <- res[[2]]
    res <- get.mutations(sons[2], edge.df, nameSeq.df, n.offspring, mutations, mutations.ind)
    mutations <- res[[1]]
    mutations.ind <- res[[2]]
  }
  return (list(mutations, mutations.ind))

}

# Calculate subtree sizes for each internal node, find mutations and compute LONR scores
#
# Arguments:
#          nameSeq.df - data frame of sequences and headers
#          edge.df - data frame with columns - parent(from), child(to), weight, distance(nt)
#          outgroup - outgroup sequence name (optional)
#
# Returns: mutations - mutation table
compute.sub.trees <- function(nameSeq.df, edge.df, outgroup = NULL){

  # remove outgroup if exists
  if (!is.null(outgroup)) {
    nameSeq.df <- nameSeq.df[-which(nameSeq.df$head == outgroup), ] # remove from sequence list
    edge.df <- edge.df[-which(edge.df[, 'to'] == outgroup), ] # remove from edge table
  }
  n.seq <- length(nameSeq.df)

  # find all subbtree roots
  roots <- as.character(unique(edge.df[!(edge.df[, 1] %in% edge.df[, 'to']), 1]))
  treesIDs <- unique(nameSeq.df$treeID)
  all.mutations <-  data.frame()
  for (i in treesIDs){
    sub.tree <- nameSeq.df[nameSeq.df$treeID==i,]
    curr.root <- as.character(roots[roots%in%sub.tree$head])

    # calculate the size of left and right subtrees for each node
    n.offspring <- data.frame(node = sub.tree$head, left = 0, right = 0, stringsAsFactors=F)
    n.offspring <- get.number.offsprings(curr.root, n.offspring, edge.df)

    # get mutations between each father-son pair
    mutations <- data.frame(mutation = rep('', 1000), LONR = 0, mutation.type = '', position=0, father='', son='', flag=F, stringsAsFactors=F)
    res <- get.mutations(curr.root, edge.df, sub.tree, n.offspring, mutations, 1)

    # remove end of table if not used
    mutations <- res[[1]]
    ind <- res[[2]]
    if (nrow(mutations) > ind )
      mutations <- mutations[-(ind:nrow(mutations)), ]

    all.mutations <- rbind(all.mutations, mutations)
  }
  return(all.mutations)

}

# Compute consensus sequence and remove positions with gaps
remove.gaps <- function(infile){

  # read FASTA file in alignment object
  aligned.seq <-read.alignment(infile, format='fasta', forceToLower = F)

  # comnpute consensus sequence
  consensus.seq <- consensus(aligned.seq, method = "majority")

  # remove columns containing gaps in consensus
  gapped.pos <- rev(which(consensus.seq == '-'))
  for (pos in gapped.pos)
    aligned.seq[['seq']] <- paste0(substr(aligned.seq[['seq']], 1, pos-1), substr(aligned.seq[['seq']], pos+1, nchar(aligned.seq[['seq']])[1]))

  # convert sequences and headers into list
  fasta.df <- data.frame(head=aligned.seq[['nam']], seq = aligned.seq[['seq']], stringsAsFactors = F)

  return (fasta.df)
}

# Dived trees into subtrees by cutting branches with more mutations than specified by cutoff
#
# Arguments:  nameSeq.df - data frame with headers and sequences
#          edge.df - data frame with columns - parent(from), child(to), edge weight (weight), edge length (distance in nt)
#          cutoff - number of mutation (default - 10)
#
# Returns: modified nameSeq.df - new column treeID, specifying sub tree
#          modified edge.df - without long branches
cut.trees <- function(nameSeq.df,edge.df,cutoff){

  # add column for tree ID
  nameSeq.df$treeID <- -1

  # Cut edges longer than threshold
  trim.edge.df <- edge.df[edge.df$distance<=cutoff,]

  # Get roots of all subtrees
  sub.roots <- setdiff(trim.edge.df$from, trim.edge.df$to)
  tree.id = 1
  # For each root, get subtree edges and nodes
  trees <- list()
  for(root in sub.roots) {
    at.bottom <- F
    nodes <- nameSeq.df[nameSeq.df$head==root,]
    # Trace down tree until at all leaves
    while(at.bottom != T) {
      # Add children of all nodes in tree thus far
      new.nodes <- unique(rbind(nodes, nameSeq.df[nameSeq.df$head %in% trim.edge.df$to[trim.edge.df$from %in% nodes$head],]))

      # If no children are to be added (tree is complete)
      if(nrow(nodes) == nrow(new.nodes)) {
        nodes <- new.nodes
        # delete internal sequences
        #nodes <- nodes[nodes$head  %in%  nodes$head[grep('^L', nodes$head, perl=T, fixed=F)],]
        if(nrow(nodes)!=0){
          nameSeq.df[nameSeq.df$head %in% nodes$head,'treeID'] <- tree.id
          tree.id <- tree.id + 1
        }
        at.bottom <- T
      }else{
        nodes <- new.nodes
      }
    }
  }
  sing <- which(nameSeq.df$treeID==-1)
  if (length(sing) > 0 )
    nameSeq.df <- nameSeq.df[-sing,]

  return(list(nameSeq.df, trim.edge.df))
}

# make sure the dirs have trailing slashes (all the path manipulation assumes they do)
check.dirs <- function(dirname) {
  if(nchar(dirname) == 0)
    stop('unexpected zero length dir name')

  last.char = substr(dirname, nchar(dirname), nchar(dirname) + 1)
  if(last.char != '/') {
      print(paste0('note: directory names must have trailing slashes, adding to ', dirname))
      dirname = paste0(dirname, '/')
  }

  return(dirname)
}

# MAIN function - Builds lineage tree, find mutations within the tree and compute LONR scores
#
# Arguments:
#          method - dnapars or neighbor
#          infile - input fasta file
#          workdir - temporary working directory
#          outgroup - outgroup sequence name (optional)
compute.LONR <- function(method, infile, workdir, outgroup=NULL, existing.edgefile=NULL, existing.node.seqfile=NULL, cutoff=10){
  workdir = check.dirs(workdir)

  # remove gaps in consensus
  fasta.df <- remove.gaps(infile)

  # fail on small or huge files
  n.seq <- nrow(fasta.df)
  if(n.seq < MIN.SEQ)
    stop('Not enough sequences to make tree with Phylip')
  if(n.seq > MAX.SEQ)
    stop('Too many sequences to make tree with Phylip')

  #------------------------------------
  # PART I - Build lineage tree
  #------------------------------------
  res <- build.trees(method, fasta.df, workdir, outgroup, existing.edgefile, existing.node.seqfile)

  if (is.null(res))
    return(F)
  nameSeq.df <- res[[1]]
  edge.df <- res[[2]]

  #------------------------------------
  # PART II - Compute LONR scores
  #------------------------------------
  res <- cut.trees(nameSeq.df,edge.df,cutoff)
  nameSeq.df <- res[[1]]
  edge.df <- res[[2]]


  LONR.table <- compute.sub.trees(nameSeq.df, edge.df, outgroup)

  # write lonr output to csv
  write.table(LONR.table, file=paste0(workdir, G.lonrfname), quote=F, sep=',', col.names=T, row.names=F)

  file.remove(paste0(workdir, G.phy.infname))
}
