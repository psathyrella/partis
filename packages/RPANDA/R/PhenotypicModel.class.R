#The superclass
setClass(
    Class = "PhenotypicModel",
    representation = representation(
        name= "character",
        period = "numeric",
        aAGamma = "function",
        numbersCopy = "numeric",
        numbersPaste = "numeric",
        initialCondition = "function",
        paramsNames = "character",
        constraints = "function",
        params0 = "numeric",
        tipLabels = "character",
        tipLabelsSimu = "character",
        comment = "character"
    ),
    prototype=prototype(
        name = "BMtest",
        period = c(0,1,2,3,4,5,6),
        aAGamma = function(i, params){
            functiona <- function(t){
                return(rep(0,i+1))
            }
            matrixA <- diag(0, i+1)
            functionGamma <- function(t){
                return(diag(params[1], i+1))
            }
            return(list(a=functiona, A=matrixA, Gamma=functionGamma))
        },
        numbersCopy = c(1, 1, 2, 1, 2, 5),
        numbersPaste = c(2, 3, 4, 5, 6, 7),
        initialCondition = function(params){
            return(list(mean=c(0,0), var=c(0,0)))
        },
        paramsNames = c("sigma"),
        constraints = function(params){
            return(params[1] > 0)
        },
        params0 = c(1),
        tipLabels = c("A", "B", "C", "D", "E", "F", "G"),
        tipLabelsSimu = c("A", "B", "C", "D", "E", "F", "G"),
        comment = "Toy model defined by defaut"
    ),
    validity=function(object){
        if( length(object@numbersCopy) != length(object@numbersPaste) ){
            stop("[PhenotypicModel : validation] \n The sequence of positions of branching lineages \n and the sequence of new positions for the traits in the newly born lineages\n should have the same length.")
        }
        if( length(object@numbersCopy) != length(object@period) ){
            stop("[PhenotypicModel : validation] \n The sequence of positions of branching lineages \n and the sequence of time periods \n should have the same length.")
        }
        if( length(object@params0) != length(object@paramsNames) ){
            stop("[PhenotypicModel : validation] \n There should be the same number of defaut parameters \n and parameter names.")
        }
        return(TRUE)
    }
)


###################################
#    Subclasses
###################################
#a new subclass has been defined for each class of models for which the "getTipDistribution" function had been optimized.

setClass(
    Class = "PhenotypicACDC",
    representation = representation(
        matrixCoalescenceTimes="matrix"
    ),
    contains="PhenotypicModel"
)

setClass(
    Class = "PhenotypicADiag",
    representation = representation(),
    contains="PhenotypicModel"
)

setClass(
    Class = "PhenotypicBM",
    representation = representation(
        matrixCoalescenceTimes="matrix"
    ),
    contains="PhenotypicModel"
)

setClass(
    Class = "PhenotypicDD",
    representation = representation(
        matrixCoalescenceJ="matrix",
        nLivingLineages="numeric"
    ),
    contains="PhenotypicModel"
)

setClass(
    Class = "PhenotypicGMM",
    representation = representation(
        n1="numeric",
        n2="numeric"
    ),
    contains="PhenotypicModel"
)

setClass(
    Class = "PhenotypicOU",
    representation = representation(
        matrixCoalescenceTimes="matrix"
    ),
    contains="PhenotypicModel"
)

setClass(
    Class = "PhenotypicPM",
    representation = representation(),
    contains="PhenotypicModel"
)


###################################
#    Basic methods
###################################

setMethod(
    f="[",
    signature="PhenotypicModel",
    definition=function(x,i,j,drop){
        switch( EXPR=i,
                "name"={return(x@name)},
                "period"={return(x@period)},
                "aAGamma"={return(x@aAGamma)},
                "numbersCopy"={return(x@numbersCopy)},
                "numbersPaste"={return(x@numbersPaste)},
                "initialCondition"={return(x@initialCondition)},
                "paramsNames"={return(x@paramsNames)},
                "constraints"={return(x@constraints)},
                "params0"={return(x@params0)},
                "tipLabels"={return(x@tipLabels)},
                "tipLabelsSimu"={return(x@tipLabelsSimu)},
                "comment"={return(x@comment)},
                stop("This variable name does not exist !")
        )
    }
)

setReplaceMethod(
    f="[",
    signature="PhenotypicModel",
    definition=function(x,i,j,value){
        switch( EXPR=i,
                "name"={x@name <- value},
                "period"={x@period <- value},
                "aAGamma"={x@aAGamma <- value},
                "numbersCopy"={x@numbersCopy <- value},
                "numbersPaste"={x@numbersPaste <- value},
                "initialCondition"={x@initialCondition <- value},
                "paramsNames"={x@paramsNames <- value},
                "constraints"={x@constraints <- value},
                "params0"={x@params0 <- value},
                "tipLabels"={x@tipLabels <- value},
                "tipLabelsSimu"={x@tipLabelsSimu <- value},
                "comment"={x@comment <- value},
                stop("This variable name does not exist !")
        )
        validObject(x)
        return(x)
    }
)

###################################
#    Affichage
###################################

setMethod(
    f="print",
    signature="PhenotypicModel",
    definition=function(x, ...){
        cat("****************************************************************\n")
        cat("*** Object of Class PhenotypicModel *** \n")
        cat("*** Name of the model : ")
        print(x@name)
        cat("*** Parameters of the model : ")
        print(x@paramsNames)
        cat("*** Description : ")
        cat(x@comment)
        cat(paste("\n*** Epochs : the model is cut into ", length(x@period), " parts. \n"))
        print(x@period)
        cat("*** Lineages branching (to be copied at the end of the corresponding period) :\n")
        print(x@numbersCopy)
        cat("*** Positions of the new trait at the end of each period :\n")
        print(x@numbersPaste)
        cat("*** Initial condition :\n")
        print(x@initialCondition)
        cat("*** Vectors a_i, A_i, Gamma_i on each period i : \n")
        print(x@aAGamma)
        cat("*** Constraints on the parameters : \n")
        print(x@constraints)
        cat("*** Defaut parameter values : ")
        print(x@params0)
        cat("*** Tip labels : \n")
        print(x@tipLabels)
        cat("*** Tip labels for simulations : \n")
        print(x@tipLabelsSimu)
        cat("****************************************************************\n")
    }
)

setMethod(
    f="show",
    signature="PhenotypicModel",
    definition=function(object){
        cat("****************************************************************\n")
        cat("*** Object of Class PhenotypicModel *** \n")
        cat("*** Name of the model : ")
        print(object@name)
        cat("*** Parameters of the model : ")
        print(object@paramsNames)
        cat("*** Description : ")
        cat(object@comment)
        cat(paste("\n*** Periods : the model is cut into ", length(object@period), " parts. \n"))
        cat("For more details on the model, call : print(PhenotypicModel)\n")
        cat("****************************************************************\n")
    }
)

