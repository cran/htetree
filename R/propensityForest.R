#' Propensity Forest for Causal Effect Regression and Estimation
#' (Modified Causal Tree Ensembles)
#'
#' Fit and evaluate a user selected number of \code{causalTree} models to get
#' an ensemble of \code{rpart} objects
#' Trees are split using covariates and the treatment vector instead of
#' the outcome variable, and evaluated using the complete data covariates
#' and the actual outcome variable.
#'
#' @aliases{propensityForest}
#'
#' @param formula a \link{formula}, with a response and features but
#' no interaction terms.  If this a a data frome, that is taken as
#' the model frame (see \code{\link{model.frame}}).
#'
#' @param data an optional data frame that includes the variables
#' named in the formula.
#'
#' @param treatment a vector that indicates the treatment status of each
#' observation. 1 represents treated and 0 represents control.  Only binary
#' treatment supported in this version.
#'
#' @param na.action the default action deletes all observations for which
#'   \code{y} is missing, but keeps those in which one or more predictors
#'   are missing.
#'
#' @param split.Rule causalTree splitting options, one of \code{"TOT"},
#' \code{"CT"}, \code{"fit"}, \code{"tstats"}, four splitting rules in
#' \code{causalTree}.  Note that the \code{"tstats"} alternative does not
#' have an associated cross-validation method \code{cv.option}; see Athey
#' and Imbens (2016) for a discussion.  Note further that \code{split.Rule}
#' and \code{cv.option} can mix and match.
#'
#' @param split.Honest boolean option, \code{TRUE} or \code{FALSE}, used
#' for \code{split.Rule} as \code{"CT"} or \code{"fit"}. If set as \code{TRUE},
#' do honest splitting, with default \code{split.alpha} = 0.5; if set as
#' \code{FALSE}, do adaptive splitting with \code{split.alpha} = 1.
#' The user choice of \code{split.alpha} will be ignored if \code{split.Honest}
#' is set as \code{FALSE}, but will be respected
#' if set to \code{TRUE}.  For \code{split.Rule}=\code{"TOT"}, there is no
#' honest splitting option and
#' the parameter \code{split.alpha} does not matter.  For
#' \code{split.Rule}=\code{"tstats"}, a value of \code{TRUE} enables use
#' of \code{split.alpha} in calculating the risk function, which determines
#' the order of pruning in cross-validation. Note also that causalTree function
#' returns the estimates from the training data, no matter what the
#' value of \code{split.Honest} is; the tree must be re-estimated to get
#' the honest estimates using \code{estimate.causalTree}. The wrapper
#' function \code{honest.CausalTree}
#' does honest estimation in one step and returns a tree.
#'
#' @param split.Bucket boolean option, \code{TRUE} or \code{FALSE}, used
#' to specify whether to apply the discrete method in splitting the tree.
#' If set as \code{TRUE}, in splitting a node, the observations in a leaf
#' will be be partitioned into buckets, with each bucket containing
#' \code{bucketNum} treated and \code{bucketNum} control units, and where
#' observations are ordered prior to partitioning. Splitting will take
#' place by bucket.
#'
#' @param bucketNum number of observations in each bucket when set
#' \code{split.Bucket} = \code{TRUE}.  However, the code will override
#' this choice in order to guarantee that there are at least \code{minsize}
#' and at most \code{bucketMax} buckets.
#'
#' @param bucketMax Option to choose maximum number of buckets to use in
#' splitting when set \code{split.Bucket} = \code{TRUE}, \code{bucketNum}
#' can change by choice of \code{bucketMax}.
#'
#' @param cv.option cross validation options, one of \code{"TOT"},
#' \code{"matching"}, \code{"CT"}, \code{"fit"}, four cross validation
#' methods in \pkg{causalTree}.  There is no \code{cv.option} for the
#' \code{split.Rule} \code{"tstats"}; see Athey and Imbens (2016) for
#' discussion.
#'
#' @param cv.Honest boolean option, \code{TRUE} or \code{FALSE}, only
#' used for \code{cv.option} as \code{"CT"} or \code{"fit"}, to specify
#' whether to apply honest risk evalation function in cross validation.
#' If set \code{TRUE}, use honest risk function, otherwise use adaptive
#' risk function in cross validation.  If set \code{FALSE}, the user choice
#' of \code{cv.alpha} will be set to 1.  If set \code{TRUE}, \code{cv.alpha}
#' will default to 0.5, but the user choice of \code{cv.alpha} will be
#' respected.  Note that honest cv estimates within-leaf variances and
#' may perform better with larger leaf sizes and/or small number of
#' cross-validation sets.
#'
#' @param minsize in order to split, each leaf must have at least
#' \code{minsize} treated cases and \code{minsize} control cases.
#' The default value is set as 2.
#'
#' @param propensity propensity score used in \code{"TOT"} splitting
#' and \code{"TOT"}, honest \code{"CT"} cross validation methods. The default
#' value is the proportion of treated cases in all observations.  In this
#' implementation, the propensity score is a constant for the whole
#' dataset.  Unit-specific propensity scores are not supported;
#' however, the user may use inverse propensity scores as case weights
#' if desired.
#'
#' @param control a list of options that control details of the
#'   \code{rpart} algorithm.  See \code{\link{rpart.control}}.
#'
#' @param split.alpha scale parameter between 0 and 1, used in splitting
#' risk evaluation function for \code{"CT"}. When \code{split.Honest = FALSE},
#' \code{split.alpha} will be set as 1.  For \code{split.Rule}=\code{"tstats"},
#' if \code{split.Honest}=\code{TRUE}, \code{split.alpha} is used in
#' calculating the risk function, which determines the order of
#' pruning in cross-validation.
#'
#' @param cv.alpha scale paramter between 0 and 1, used in cross
#' validation risk evaluation function for \code{"CT"} and \code{"fit"}.  When
#' \code{cv.Honest = FALSE}, \code{cv.alpha} will be set as 1.
#'
#' @param sample.size.total Sample size used to build each tree in the
#' forest (sampled randomly with replacement)
#'
#' @param sample.size.train.frac Fraction of the sample size used for
#' building each tree (training)
#'
#' @param mtry Number of data features used to build a tree (This
#' variable is not used presently)
#'
#' @param nodesize Minimum number of observations for treated and
#' control cases in one leaf node
#'
#' @param num.trees Number of trees to be built in the causal forest
#'
#' @param ncolx Total number of covariates
#'
#' @param ncov_sample Number of covariates randomly sampled to build
#' each tree in the forest
#'
#'
#' @details
#' Propensity forest is similar to a causal forest, with some important
#' differences as discussed below. The causalForest builds an ensemble of
#' CausalTrees, by repeated random sampling of the data with replacement.
#' For prediction, the average value over all tree predictions is used.
#' Propensity forest differs from a causal forest in the following way:
#' The tree building phase is done by using the covariates and treatment
#' vector as the dummy output/outcomes variable (instead of the actual
#' outcomes variable) During the tree evaluation phase, the
#' reestimation error is calculated on the actual outcomes variable to
#' evaluate the tree performance.
#' Note that the propensity forest will always build an adaptive
#' (non honest) ensemble of trees.
#' CausalTree differs from \code{rpart} function from \pkg{rpart}
#' package in splitting rules and cross validation methods. Please check
#' Athey and Imbens, \emph{Recursive Partitioning for Heterogeneous Causal
#' Effects} (2016) and Stefan Wager and Susan Athey, \emph{Estimation
#' and Inference of Heterogeneous Treatment Effects using Random Forests
#' } for more details.
#' @returns An object of class \code{rpart}.  See \code{\link{rpart.object}}.
#'
#' @examples
#' library(rpart)
#' library("htetree")
#' pf <- propensityForest(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simulation.1,
#'   treatment=simulation.1$treatment,
#'   split.Bucket=FALSE,
#'   sample.size.total = floor(nrow(simulation.1) / 2),
#'   sample.size.train.frac = .5,
#'   mtry = ceiling(ncol(simulation.1)/3), nodesize = 25, num.trees= 5,
#'   ncolx=10,ncov_sample=3)
#'
#' pfpredtest <- predict(pf, newdata=simulation.1[1:100,],
#'   type="vector")

propensityForest <- function(formula, data, treatment,
     na.action = na.causalTree,
     split.Rule="CT", split.Honest=TRUE, split.Bucket=FALSE, bucketNum = 5,
     bucketMax = 100, cv.option="CT", cv.Honest=TRUE, minsize = 2L,
     propensity=mean(treatment), control, split.alpha = 0.5, cv.alpha = 0.5,
     sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = 1,
     mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),
     ncolx, ncov_sample) {

  # do not implement subset option of causalTree, inherited from rpart
  # do not implement weights and costs yet

  if(sample.size.train.frac != 1) {
    # print("warning: for propensity Forest, sample.size.train.frac should be 1; resetting to 1")
    sample.size.train.frac <- 1
  }

  num.obs <-nrow(data)

  vars <- all.vars(formula)
  y <- vars[1]
  x <- vars[2:length(vars)]
  treatmentdf <- data.frame(treatment)
  # if(class(data)[1]=="data.table"){
  if(inherits(data,"data.table",TRUE)==1){
    treatmentdt <- data.table(treatment)
    datax<-data[,..x]
    datay<-data[,y,with=FALSE]
    data <- cbind(datax,datay, treatmentdt)
  # }else if(class(data)=="data.frame"){
  }else if(inherits(data,"data.frame")){
    data <- data[, c(x, y)]
    data <- cbind(data, treatmentdf)
  }

  causalForest.hon <- init.causalForest(formula=formula, data=data, treatment=treatment, num.trees=num.trees, weights=FALSE, cost=FALSE,ncov_sample=ncov_sample)
  sample.size <- min(sample.size.total, num.obs)
  train.size <- round(sample.size.train.frac*sample.size)

  outcomename = as.character(formula[2])

  # print("Building trees ...")

  for (tree.index in 1:num.trees) {

    # print(paste("Tree", as.character(tree.index)))

    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    train.idx <- full.idx[1:train.size]

    cov_sample<-sample.int(ncolx)
    cov_sample<-cov_sample[1:ncov_sample]

    #modify the y=f(x) equation accordingly for this tree
    #and modify the colnames
    fsample<-""
    nextx<-""
    nameall_sample<-c()
    for (ii in 1:(ncov_sample)){
      nextxindex <- cov_sample[ii]
      nextx <- x[[nextxindex]]
      if (ii==1) {
        fsample <-nextx
        name<- nextx
      }
      if (ii>1) {
        fsample <- paste0(fsample,"+",nextx)
        name <- c(name, nextx)
      }
    }

    #nameall_sample <- c( name,"temptemp","y", "tau_true","treattreat")
    nameall_sample <- c( name,"temptemp",y, "treattreat")
    nameall_sample_save <- c( name,  y, "w") #, "tau_true")

    #store this var subset for each tree (need it during testing/predict stage)
    causalForest.hon$cov_sample[tree.index,]<-cov_sample
    #also store the formula & colnames of X for each tree (need it during testing/predict stage)
    causalForest.hon$nameall_sample[tree.index,]<-nameall_sample_save
    causalForest.hon$fsample[[tree.index]]<-fsample

    # rename variables as a way to trick rpart into building the tree with all the object attributes considering the outcome variable as named
    # by the input formula, even though the tree itself is trained on w.  Note that we aren't saving out this propensity tree anyway, but if
    # we decided later to try to save out the propensity trees and do something directly with the propensity scores, we would need to do something
    # more tedious like estimate the propensity tree with the original names, and then edit the attributes to replace the treatment variable name
    # with the outcome variable name for the estimate part
    dataTree <- data.frame(data[train.idx,])
    dataTree$treattreat <- treatmentdf[train.idx,]
    names(dataTree)[names(dataTree)==outcomename] <- "temptemp"
    names(dataTree)[names(dataTree)=="treattreat"] <- outcomename


    # if(class(data)[1]=="data.table"){
    if(inherits(data,"data.table",TRUE)==1){
      #sample covariates and pick relevant covariates for tree
      treeRange<-c(cov_sample,(ncolx+1):ncol(dataTree))
      dataTree <- dataTree[,..treeRange]
    # }else if(class(data)=="data.frame"){
    }else if(inherits(data,"data.frame")){
      #pick relevant covariates for tree
      dataTree <- dataTree[,c(cov_sample,(ncolx+1):ncol(dataTree))]
    }

    #change colnames to reflect the sampled cols
    names(dataTree)=nameall_sample
    # names(dataEstim)=nameall_sample
    formula<-paste(y,"~",fsample,sep="")
    #one options: estimate the propensity tree with anova so that it will be type "anova" when we re-estimate
    #here: replace elements of the rpart object to make it look like anova tree, so that we'll be able to properly predict with it later, etc.
    tree.propensity <- rpart(formula=formula, data=dataTree, method="class",
                             control=rpart.control(cp=0, minbucket=nodesize))

    # make it look like a method="anova" tree
    tree.propensity$method <- "anova"
    tree.propensity$frame$yval2 <- NULL
    tree.propensity$functions$print <- NULL

    # switch the names back in the data frame so that when we estimate treatment effects, will have the right outcome variables
    names(dataTree)[names(dataTree)==y] <- "treattreat"
    names(dataTree)[names(dataTree)=="temptemp"] <- y
    tree.treatment <- estimate.causalTree(object=tree.propensity,data=dataTree, treatment=dataTree$treattreat)

    causalForest.hon$trees[[tree.index]] <- tree.treatment
    causalForest.hon$inbag[full.idx, tree.index] <- 1
  }

  return(causalForest.hon)
}
