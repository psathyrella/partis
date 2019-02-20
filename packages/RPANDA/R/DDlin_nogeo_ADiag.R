##################################################
#    Bank of Classic 1D Phenotypic Models
##################################################

.createModel_DDlin <- function(tree){
    

        comment <- "Diversity dependent exponential model \n Implemented as in Drury et al. Systematic Biology."
        paramsNames <- c("m0", "logsigma0", "r")
        params0 <- c(0,log(1),-0.1)

		###probably needs to be updated to reflect times/spans in geo object, perhaps could feed directly to model call below

        periodizing <- periodizeOneTree(tree) 
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        
        ##need to figure out what this does, should I force it to 0 or do I need to call VCV matrix to get GLS solution for mean?
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
        
            
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(rep(0, length(vectorU)))
            matrixGamma <- function(t) return(sqrt((exp(params[2])^2)+(params[3]*length(vectorU)))*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return((exp(params[2])^2)+(params[3]*length(tree$tip.label)) > 0)

        model <- new(Class="PhenotypicADiag", name="DDlin", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling, comment=comment)


    return(model)
}



##################################################
#    Describe the periods on a 'phylo' tree
##################################################

getMatrixCoalescenceJ <- function(tree, periods){
    # The entry (k,l) of the matrix is the index j such that tau_j = t_{k,l}
    matrixCoalescenceTimes <- findMRCA(tree, type="height")
    n <- length(matrixCoalescenceTimes[,1])
    matrixCoalescenceJ <- diag(0, n)
    for(k in 1:n){
        for(l in 1:n){
            matrixCoalescenceJ[k,l] <- which(periods == matrixCoalescenceTimes[k,l])
        }
    }

    return(matrixCoalescenceJ)
}

isATip <- function(tree, branch_number){
    return(!(tree$edge[branch_number,2] %in% tree$edge[,1]))
}

periodizeOneTree <- function(tree){
    # Returns 3 vectors giving 
    # 1) the periods of the tree, 
    # 2) the starting times of all branches in the tree 
    # 3) the death time of all branches in the tree
    
    nodeheight <- nodeHeights(tree)
    startingTimes <- nodeheight[,1]
    endTimes <- nodeheight[,2]
    all_time_events <- sort(c(startingTimes, endTimes))
    # the following removes identical entries in the vector
    periods <- unique(all_time_events)
    
    return(list(periods=periods, startingTimes=startingTimes, endTimes=endTimes))
}

endOfPeriods <- function(periodizing, tree){
    # Returns the list of branching or dying lineages at the beginning of each period : copy
    # Together with the list of places where the new lineage is inserted (or zero if a lineage dies) : paste
    # And the number of lineages on the focal period : nLineages
    # The rule is : at each branching point, the first of the two new branches is assigned its mother label, and the new one takes the last label (n, where n is the number of lineages at that time)
    
    nBranch <- length(periodizing$startingTimes)
    nPeriods <- length(periodizing$periods)
    
    numbersCopy <- rep(0, times=nPeriods)
    numbersPaste <- rep(0, times=nPeriods)
    numbersLineages <- rep(0, times=nPeriods)
    numbersLivingLineages <- rep(0, times=nPeriods)

    # We initialize the labeling of branches in the tree
    labelingLineages <- rep(0, times=nBranch)
    initialBranches <- periodizing$startingTimes[periodizing$startingTimes==0]
    if(length(initialBranches) == 1){
        labelingLineages[1] <- 1
        n <- 1
    }else{
        labelingLineages[periodizing$startingTimes==0] <- c(1,2)
        n <- 2
    }
    numbersLineages[1] <- n
    numbersLivingLineages[1] <- n
    numbersCopy[1] <- 1
    numbersPaste[1] <- 2
    
    for(i in 2:nPeriods){
        tau_i <- periodizing$periods[i]
        newBranches <- which(tau_i == periodizing$startingTimes)
        # If tau_i is a birth time on the tree
        if(length(newBranches) == 2){
            n <- n+1
            labelingLineages[newBranches[1]] <- labelingLineages[newBranches[1]-1]
            labelingLineages[newBranches[2]] <- n
            numbersCopy[i] <- labelingLineages[newBranches[1]-1]
            numbersPaste[i] <- n
            numbersLivingLineages[i] <- numbersLivingLineages[i-1]+1
        # Else, tau_i is only a death time of one or many terminal branches.
        }else{
            deadBranches <- which(tau_i == periodizing$endTimes)
            numbersCopy[i] <- labelingLineages[ deadBranches[1] ]
            numbersPaste[i] <- 0
            numbersLivingLineages[i] <- numbersLivingLineages[i-1]-1
        }
        numbersLineages[i] <- n
    }

    permutationLabels <- labelingLineages[!(periodizing$endTimes %in% periodizing$startingTimes)]
    labeling <- tree$tip.label[order(permutationLabels)]
    
    return(list(copy=numbersCopy, paste=numbersPaste, nLineages=numbersLineages, labeling=labeling, nLivingLineages=numbersLivingLineages))
}

getLivingLineages <- function(i, eventEndOfPeriods){
    
    livingLineages <- rep(1, times=eventEndOfPeriods$nLineages[i])
    deads <- eventEndOfPeriods$copy[1:i][eventEndOfPeriods$paste[1:i] == 0]
    livingLineages[deads] <- 0
    
    return(livingLineages)
}