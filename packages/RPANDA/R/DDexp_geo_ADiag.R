.createModel_DDexp_geo <- function(tree,geo.object){
    

        comment <- "Diversity dependent exponential model with biogeography\n Implemented as in Drury et al. Systematic Biology."
        paramsNames <- c("m0", "logsigma0", "r")
        params0 <- c(0,log(1),-0.1)

        periodizing <- periodizeOneTree_geo(tree,geo.object) 
        eventEndOfPeriods <- endOfPeriods(periodizing, tree)
        
        initialCondition <- function(params) return( list(mean=c(params[1]), var=matrix(c(0))) ) 
        
            
        aAGamma <- function(i, params){
            vectorU <- getLivingLineages(i, eventEndOfPeriods)
            vectorA <- function(t) return(rep(0, length(vectorU)))
            matrixGamma <- function(t) return(exp(params[2])*exp((params[3]/2)*colSums(geo.object$geography.object[[i]]))*diag(vectorU))
            matrixA <- diag(0, length(vectorU))
            return(list(a=vectorA, A=matrixA, Gamma=matrixGamma))
        }

        constraints <- function(params) return(params[3]<=Inf && params[3] >= -Inf)
        
        model <- new(Class="PhenotypicADiag", name="DDexp_geo", period=periodizing$periods, aAGamma=aAGamma, numbersCopy=eventEndOfPeriods$copy, numbersPaste=eventEndOfPeriods$paste, initialCondition=initialCondition, paramsNames=paramsNames, constraints=constraints, params0=params0, tipLabels=eventEndOfPeriods$labeling,  comment=comment)

    return(model)
}



##################################################
#    Describe the periods on a 'phylo' tree
##################################################

getStartingTimes <- function(tree){
    # Returns a vector giving the starting time for each branch of a tree in format "phylo"
    
    nBranch = length(tree$edge.length)
    starting_times <- rep(0, times=nBranch)
    
    # we add progressively for each branch the length of all parent branches in the vector "starting_times"
    for(n1 in 1:nBranch){
        n2 <- n1 + 1
        while(n2 <= nBranch){
            if(tree$edge[n2,1]==tree$edge[n1,2]){
                starting_times[n2] <- starting_times[n1] + tree$edge.length[n1]
            }
            n2 <- n2+1
        }
    }
    
    return(starting_times)
}


isATip <- function(tree, branch_number){
    return(!(tree$edge[branch_number,2] %in% tree$edge[,1]))
}

periodizeOneTree_geo <- function(tree,geo.object){
    # Returns 3 vectors giving 
    # 1) the periods of the tree, 
    # 2) the starting times of all branches in the tree 
    # 3) the death time of all branches in the tree
    hold<-nodeHeights(tree)
    startingTimes <- hold[,1]
    endTimes <- hold[,2]
    all_time_events <- sort(c(startingTimes, endTimes))
    
    
    nodetimes=max(branching.times(tree))-sort(branching.times(tree),decreasing=TRUE)
	extv<-vapply(geo.object$geography.object,function(x)dim(x)[1],1)
	outv<-c(1)
	for(i in 2:length(extv)){
		if(extv[i]!=extv[i-1]){
			outv<-c(outv,i)
		}}
	
	chg.times=which(!1:length(geo.object$times)%in%c(outv,length(geo.object$times)))
	periods=sort(c(geo.object$times[chg.times],unique(startingTimes),max(endTimes)))
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
        # Else, tau_i is only a death time of one or many terminal branches.
        }else{
            numbersCopy[i] <- 0
 			#deadBranches <- which(tau_i == periodizing$endTimes)
            #numbersCopy[i] <- labelingLineages[ deadBranches[1] ]
            numbersPaste[i] <- 0
        }
        numbersLineages[i] <- n
    }

    permutationLabels <- labelingLineages[!(periodizing$endTimes %in% periodizing$startingTimes)]
    labeling <- tree$tip.label[order(permutationLabels)]
    
    return(list(copy=numbersCopy, paste=numbersPaste, nLineages=numbersLineages, labeling=labeling))
}


getLivingLineages <- function(i, eventEndOfPeriods){
    
    livingLineages <- rep(1, times=eventEndOfPeriods$nLineages[i])
    deads <- eventEndOfPeriods$copy[1:i][eventEndOfPeriods$paste[1:i] == 0]
    livingLineages[deads] <- 0
    
    return(livingLineages)
}
    