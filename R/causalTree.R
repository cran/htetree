#' Causal Effect Regression and Estimation Trees
#'
#' Fit a \code{causalTree} model to get an \code{rpart} object
#'
#' @param formula a \link{formula}, with a response and features but
#' no interaction terms.  If this a a data frome, that is taken as
#' the model frame (see \code{\link{model.frame}).}
#' @param data an optional data frame that includes the variables
#' named in the formula.
#' @param weights optional case weights.
#' @param treatment a vector that indicates the treatment status of
#' each observation. 1 represents treated and 0 represents control.
#' Only binary treatment supported in this version.
#' @param subset optional expression saying that only a subset of the
#' rows of the data should be used in the fit.
#' @param na.action the default action deletes all observations for which
#' \code{y} is missing, but keeps those in which one or more predictors
#' are missing.
#' @param split.Rule causalTree splitting options, one of \code{"TOT"},
#' \code{"CT"}, \code{"fit"}, \code{"tstats"}, four splitting rules in
#' \code{causalTree}.  Note that the \code{"tstats"} alternative does
#' not have an associated cross-validation method \code{cv.option};
#' see Athey and Imbens (2016) for a discussion.  Note further
#' that \code{split.Rule} and \code{cv.option} can mix and match.
#'
#' @param split.Honest boolean option, \code{TRUE} or \code{FALSE},
#' used for \code{split.Rule} as \code{"CT"} or \code{"fit"}.
#' If set as \code{TRUE}, do honest splitting, with default
#' \code{split.alpha} = 0.5; if set as \code{FALSE}, do adaptive
#' splitting with \code{split.alpha} = 1.  The user choice of
#' \code{split.alpha} will be ignored if \code{split.Honest} is set
#' as \code{FALSE}, but will be respected
#' if set to \code{TRUE}.  For \code{split.Rule}=\code{"TOT"},
#' there is no honest splitting option and
#' the parameter \code{split.alpha} does not matter.
#' For \code{split.Rule}=\code{"tstats"}, a value of \code{TRUE}
#' enables use of \code{split.alpha} in calculating the risk function,
#' which determines the order of pruning in cross-validation. Note also
#' that causalTree function returns the estimates from the training data,
#' no matter what the value of \code{split.Honest} is; the tree must be
#' re-estimated to get the honest estimates using \code{estimate.causalTree}.
#' The wrapper function \code{honest.CausalTree}
#' does honest estimation in one step and returns a tree.
#'
#' @param HonestSampleSize number of observations anticipated to be used
#' in honest re-estimation after building the tree. This enters the
#' risk function used in both splitting and cross-validation.
#'
#' @param split.Bucket boolean option, \code{TRUE} or \code{FALSE},
#' used to specify whether to apply the discrete method in splitting the tree.
#' If set as \code{TRUE}, in splitting a node, the observations in a leaf will
#' be be partitioned into buckets, with each bucket containing \code{bucketNum}
#' treated and \code{bucketNum} control units, and where observations are
#' ordered prior to partitioning. Splitting will take place by bucket.
#'
#' @param bucketNum number of observations in each bucket when set
#' \code{split.Bucket} = \code{TRUE}.  However, the code will override
#' this choice in order to guarantee that there are at least \code{minsize}
#' and at most \code{bucketMax} buckets.
#'
#' @param bucketMax Option to choose maximum number of buckets to use in
#' splitting when set \code{split.Bucket} = \code{TRUE},
#' \code{bucketNum} can change by choice of \code{bucketMax}.
#'
#' @param cv.option cross validation options, one of \code{"TOT"},
#' \code{"matching"}, \code{"CT"}, \code{"fit"}, four cross validation
#' methods in \pkg{causalTree}.  There is no \code{cv.option} for the
#' \code{split.Rule} \code{"tstats"}; see Athey and Imbens (2016) for
#' discussion.
#'
#' @param cv.Honest boolean option, \code{TRUE} or \code{FALSE},
#' only used for \code{cv.option} as \code{"CT"} or \code{"fit"},
#' to specify whether to apply honest risk evalation function in cross
#' validation. If set \code{TRUE}, use honest risk function, otherwise use
#' adaptive risk function in cross validation.  If set \code{FALSE}, the user
#' choice of \code{cv.alpha} will be set to 1.  If set \code{TRUE},
#' \code{cv.alpha} will default to 0.5, but the user choice of
#' \code{cv.alpha} will be respected.  Note that honest cv estimates
#' within-leaf variances and may perform better with larger leaf sizes
#' and/or small number of cross-validation sets.
#'
#' @param minsize in order to split, each leaf must have at least
#' \code{minsize} treated cases and \code{minsize} control cases.
#' The default value is set as 2.
#'
#' @param x keep a copy of the \code{x} matrix in the result.
#'
#' @param y keep a copy of the dependent variable in the result.  If
#' missing and \code{model} is supplied this defaults to \code{FALSE}.
#'
#' @param propensity propensity score used in \code{"TOT"} splitting
#' and \code{"TOT"}, honest \code{"CT"} cross validation methods. The
#' default value is the proportion of treated cases in all observations.
#' In this implementation, the propensity score is a constant for the whole
#'  dataset.  Unit-specific propensity scores are not supported; however,
#' the user may use inverse propensity scores as case weights if desired.
#'
#' @param control a list of options that control details of the
#' \code{rpart} algorithm.  See \code{\link{rpart.control}}.
#'
#' @param split.alpha scale parameter between 0 and 1, used in splitting
#' risk evaluation function for \code{"CT"}. When \code{split.Honest = FALSE},
#' \code{split.alpha} will be set as 1.  For \code{split.Rule}=\code{"tstats"},
#' if \code{split.Honest}=\code{TRUE}, \code{split.alpha} is used in
#' calculating the risk function, which determines the order of pruning
#' in cross-validation.
#'
#' @param cv.alpha scale paramter between 0 and 1, used in cross validation
#' risk evaluation function for \code{"CT"} and \code{"fit"}.  When
#' \code{cv.Honest = FALSE}, \code{cv.alpha} will be set as 1.
#'
#' @param cost a vector of non-negative costs, one for each variable in
#' the model. Defaults to one for all variables. These are scalings to
#' be applied when considering splits, so the improvement on splitting
#' on a variable is divided by its cost in deciding which split to
#' choose.
#'
#' @param cv.gamma,split.gamma optional parameters used in evaluating policies.
#'
#' @param \dots arguments to \code{\link{rpart.control}} may also be
#' specified in the call to \code{causalTree}.  They are checked against the
#' list of valid arguments.  An example of a commonly set parameter would
#' be \code{xval}, which sets the number of cross-validation samples.
#' The parameter \code{minsize} is implemented differently in
#' \code{causalTree} than in \code{rpart}; we require a minimum of \code{minsize}
#' treated observations and a minimum of \code{minsize} control
#' observations in each leaf.
#'
#' @details CausalTree differs from \code{rpart} function from \pkg{rpart}
#' package in splitting rules and cross validation methods. Please check
#' Athey and Imbens, \emph{Recursive Partitioning for Heterogeneous Causal
#' Effects} (2016) for more details.
#' @returns An object of class \code{rpart}.  See \code{\link{rpart.object}}.
#' @references Breiman L., Friedman J. H., Olshen R. A., and Stone, C. J. (1984)
#' \emph{Classification and Regression Trees.}
#' Wadsworth.
#'
#' Athey, S and G Imbens (2016)  \emph{Recursive Partitioning for Heterogeneous
#' Causal Effects}.  http://arxiv.org/abs/1504.01132
#'
#' @seealso \code{\link{honest.causalTree}},
#' \code{\link{rpart.control}}, \code{\link{rpart.object}},
#' \code{\link{summary.rpart}}, \code{\link{rpart.plot}}
#'
#' @examples
#' library("htetree")
#' library("rpart")
#' library("rpart.plot")
#' tree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1,
#' treatment = simulation.1$treatment,
#' split.Rule = "CT", cv.option = "CT", split.Honest = TRUE, cv.Honest = TRUE,
#' split.Bucket = FALSE, xval = 5,
#' cp = 0, minsize = 20, propensity = 0.5)
#'
#' opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
#'
#' opfit <- prune(tree, opcp)
#'
#' rpart.plot(opfit)
#'
#' fittree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1,
#'                       treatment = simulation.1$treatment,
#'                       split.Rule = "fit", cv.option = "fit",
#'                       split.Honest = TRUE, cv.Honest = TRUE, split.Bucket = TRUE,
#'                       bucketNum = 5,
#'                       bucketMax = 200, xval = 10,
#'                       cp = 0, minsize = 20, propensity = 0.5)
#'
#' tstatstree <- causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1,
#'                          treatment = simulation.1$treatment,
#'                          split.Rule = "tstats", cv.option = "CT",
#'                          cv.Honest = TRUE, split.Bucket = TRUE,
#'                          bucketNum = 10,
#'                          bucketMax = 200, xval = 5,
#'                          cp = 0, minsize = 20, propensity = 0.5)




causalTree <- function(formula, data, weights, treatment, subset,
					   na.action = na.causalTree,
					   split.Rule, split.Honest, HonestSampleSize, split.Bucket, bucketNum = 5,
					   bucketMax = 100, cv.option, cv.Honest, minsize = 2L,
					   x = FALSE, y = TRUE, propensity, control, split.alpha = 0.5,
					   cv.alpha = 0.5,cv.gamma=0.5,split.gamma=0.5,
					   cost, ...){
  # model.response<-model.weights <-model.offset <- .getXlevels <- rpart.control <-
    # mlist <- NULL

	Call <- match.call()

	indx <- match(c("formula", "data", "weights", "subset"),
				  names(Call), nomatch = 0L)

	if (indx[1] == 0L) stop("a 'formula' argument is required")
	temp <- Call[c(1L, indx)]
	temp$na.action <- na.action
	temp[[1L]] <- quote(stats::model.frame)
	names(treatment) <- rownames(data)
	m <- eval.parent(temp)
	treatment <- treatment[(rownames(m))]

	Terms <- attr(m, "terms")
	if (any(attr(Terms, "order") > 1L))
		stop("Trees cannot handle interaction terms")

	Y <- model.response(m)
	wt <- model.weights(m)
	if (any(wt < 0)) stop("negative weights not allowed")
	if (!length(wt)) wt <- rep(1, nrow(m))
	offset <- model.offset(m)
	X <- causalTree.matrix(m)
	#save("X",file="test.Rdata")
	nobs <- nrow(X)
	nvar <- ncol(X)
	#treatment <- m$`(treatment)`

	## requirement for treatment status
	if (missing(treatment)) {
		stop("You should input the treatment status vector for data:
			 1 represent treated and 0 represent controlled.")
	}
	if (sum(treatment %in% c(0,1)) != nobs) {
		stop("The treatment status should be 1 or 0 only: 1 represent treated and 0 represent controlled.")
	}

	if (sum(treatment) == 0 || sum(treatment) == nobs) {
		stop("The data only contains treated cases or controlled cases, please check 'treatment' again.")
	}

    ## check propensity score
	if (missing(propensity)) {
		propensity <- sum(treatment) / nobs
	}

	## check the Split.Rule:
	if (missing(split.Rule)) {
		split.Rule <- "TOT"
		warning("The default split rule is 'TOT'.")
	}

	# check split.Bucket:
	if (missing(split.Bucket)) {
		split.Bucket <- FALSE
		bucketNum <- 0
		bucketMax <- 0
	}
	split.Bucket.num <- pmatch(split.Bucket, c(TRUE, FALSE))
	if (is.na(split.Bucket.num))
		stop("Invalid split.Bucket input, split.Bucket can be only TRUE or FALSE.")

	if (!split.Bucket) {
		# split.Bucket = F
		bucketNum <- 0
		bucketMax <- 0
	} else {
		# split.Bucket = T
		if (missing(bucketMax)) {
			# maximum number of buckets
			bucketMax <- 100
		}
		if (missing(bucketNum)) {
			# number of treat cases or control cases in one bucket
			#Numbuckets = max(minsize, min(round(numtreated/bucketNum),round(numcontrol/bucketNum),bucketMax))
			bucketNum <- 5
			# at least 5 obs in one bucket
		}
		split.Rule <- paste(split.Rule, 'D', sep = '')
	}

	split.Rule.int <- pmatch(split.Rule, c("TOT", "CT", "fit", "tstats", "TOTD", "CTD",
	                                       "fitD", "tstatsD", "user", "userD","policy","policyD"))
	# print(split.Rule.int)
	# print(split.Rule)
	if (is.na(split.Rule.int)) stop("Invalid splitting rule.")
	split.Rule <- c("TOT", "CT", "fit", "tstats", "TOTD", "CTD", "fitD",
	                "tstatsD", "user", "userD","policy","policyD")[split.Rule.int]

	## check the Split.Honest, for convenience
	if (split.Rule.int %in% c(1, 5)) {
		if (!missing(split.Honest)) {
			warning("split.Honest is not used in your chosen splitting rule.")
		}
		if (!missing(split.alpha)) {
			warning("split.alpha is not used in your chosen splitting rule. split.Honest set to FALSE")
		}
		split.Honest <- FALSE
	} else {
		if (missing(split.Honest)) {
			split.Honest <- TRUE
			warning("The default split.Honest = TRUE for your chosen splitting rule.")
		}

	}

	## check the Split.Honest == T/F
	split.Honest.num <- pmatch(split.Honest, c(T, F))
	if(is.na(split.Honest.num))
		stop("Invalid split.Honest input, split.Honest can be only TRUE or FALSE.")

	if (split.Honest == TRUE && split.Rule.int %in% c(2, 3, 4, 6, 7, 8, 9, 10,11,12)) {
		# ct, fit, tstats, ctd, fitd, tstatsd, user, userd,policy,policyD:
		if(missing(split.alpha)) {
			# set default honest splitting alpha to 0.5
			split.alpha <- 0.5
		} else {
			# check split.alpha in [0, 1]
			if (split.alpha > 1 || split.alpha < 0) {
				stop("Invalid input for split.alpha. split.alpha should between 0 and 1.")
			}
		}
	  #check for gamma for policy
	  if(missing(split.gamma)) {
	    # set default honest splitting alpha to 0.5
	    split.gamma <- 0.5
	  } else {
	    # check split.alpha in [0, 1]
	    if (split.gamma > 1 || split.gamma < 0) {
	      stop("Invalid input for split.gamma. split.gamma should between 0 and 1.")
	    }
	  }

	} else if (split.Rule.int %in% c(2, 3, 4, 6, 7, 8, 9, 10,11,12)){
		# split.Honest = False
		if (split.alpha != 1)
			warning("For dishonest(adaptive) splitting, split.alpha =  1.");
		split.alpha <- 1
		#added gamma check for policy
		if(missing(split.gamma)) {
		  # set default honest splitting alpha to 0.5
		  split.gamma <- 0.5
		} else {
		  # check split.alpha in [0, 1]
		  if (split.gamma > 1 || split.gamma < 0) {
		    stop("Invalid input for split.gamma. split.gamma should between 0 and 1.")
		  }
		}
	}


	# check propensity score:
	if (split.Rule %in% c("TOT", "TOTD", "fit", "fitD")) {
		if (propensity > 1 || propensity < 0) {
			stop("Propensity score should be between 0 and 1.")
		}
	}

	xvar <- apply(X, 2, var)
	method <- "anova"
	method.int <- 1

	# ------------------------------------- cv begins -------------------------------------- #
	if (missing(cv.option)) {
		# temporarily, no crossvalidation
		warning("Miss 'cv.option', choose not to do cross validations.")
		cv.option <-"none"
		xval <- 0
	}

	if(missing(cv.Honest)) {
		cv.Honest <- TRUE
	}
	cv.Honest.num <- pmatch(cv.Honest, c(T, F))
	if (is.na(cv.Honest.num))
		stop ("Invalid cv.Honest. cv.Honest should be TRUE or FALSE.")

	if (cv.option == "CT" || cv.option == "fit" || cv.option == "user" || cv.option == "policy") {
		if (cv.Honest) {
			# cv.Honest = T
			cv.option <- paste(cv.option, 'H', sep = '')
		} else {
			cv.option <- paste(cv.option, 'A', sep = '')
		}
	}

	cv.option.num <- pmatch(cv.option, c("TOT", "matching", "fitH", "fitA", "CTH", "CTA", "userH", "userA","policyH","policyA", "none"))
	if(is.na(cv.option.num)) stop("Invalid cv option.")

	# check cv.alpha
	if (cv.option.num %in% c(1, 2, 4, 6, 8)) {
	    # tot, matching, fitA, CTA, userA
		if (!missing(cv.alpha))
			warning("cv.alpha is not used in your chosen cross validation method.")
	}

	if (missing(cv.alpha)) {
		cv.alpha <- 0.5
	}

	#for policy, set gamma (set for all presently)
	if (missing(cv.gamma)) {
	  cv.gamma <- 0.5
	}
	if (missing(HonestSampleSize)) {
		# default should be the # of samples in training
		HonestSampleSize <- nobs
	}

	HonestSampleSize <- as.integer(HonestSampleSize)
	if (is.na(HonestSampleSize))
		stop("HonestSampleSize should be an integer.")
	# -------------------------------- cv checking ends -------------------------------- #

	init <- get(paste("htetree", method, sep = "."), envir = environment())(Y, offset, wt)

	ns <- asNamespace("htetree")
	if (!is.null(init$print)) environment(init$print) <- ns
	if (!is.null(init$summary)) environment(init$summary) <- ns
	if (!is.null(init$text)) environment(init$text) <- ns

	Y <- init$y

	xlevels <- .getXlevels(Terms, m)
	cats <- rep(0L, ncol(X))
	if (!is.null(xlevels))
		cats[match(names(xlevels), colnames(X))] <-
			unlist(lapply(xlevels, length))

		extraArgs <- list(...)

		if (length(extraArgs)) {
			controlargs <- names(formals(rpart.control)) # legal arg names
			indx <- match(names(extraArgs), controlargs, nomatch = 0L)

			if (any(indx == 0L))
				stop(gettextf("Argument %s not matched",
							  names(extraArgs)[indx == 0L]),
					 domain = NA)
		}

		controls <- causalTree.control(...)
		if (!missing(control)) controls[names(control)] <- control

		xval <- controls$xval
		if (is.null(xval) || (length(xval) == 1L && xval == 0L) || method=="user") {
			xgroups <- 0L
			xval <- 0L
		} else if (length(xval) == 1L) {
			## make random groups
			control_idx <- which(treatment == 0)
			treat_idx <- which(treatment == 1)
			xgroups <- rep(0, nobs)
			xgroups[control_idx] <- sample(rep(1L:xval, length = length(control_idx)), length(control_idx), replace = F)
			xgroups[treat_idx] <- sample(rep(1L:xval, length = length(treat_idx)), length(treat_idx), replace = F)
		} else if (length(xval) == nobs) {
			## pass xgroups by xval
			xgroups <- xval
			xval <- length(unique(xgroups))
		} else {
			## Check to see if observations were removed due to missing
			if (!is.null(attr(m, "na.action"))) {
				## if na.causalTree was used, then na.action will be a vector
				temp <- as.integer(attr(m, "na.action"))
				xval <- xval[-temp]
				if (length(xval) == nobs) {
					xgroups <- xval
					xval <- length(unique(xgroups))
				} else stop("Wrong length for 'xval'")
			} else stop("Wrong length for 'xval'")
		}


		##
		## Incorprate costs
		##
		if (missing(cost)) cost <- rep(1, nvar)
		else {
			if (length(cost) != nvar)
				stop("Cost vector is the wrong length")
		if (any(cost <= 0)) stop("Cost vector must be positive")
		}

		##
		## Have C code consider ordered categories as continuous
		##  A right-hand side variable that is a matrix forms a special case
		## for the code.
		##
		tfun <- function(x)
			if (is.matrix(x)) rep(is.ordered(x), ncol(x)) else is.ordered(x)

		labs <- sub("^`(.*)`$", "\\1", attr(Terms, "term.labels")) # beware backticks
		isord <- unlist(lapply(m[labs], tfun))

		storage.mode(X) <- "double"
		storage.mode(wt) <- "double"
		storage.mode(treatment) <- "double"
		minsize <- as.integer(minsize) # minimum number of obs for treated and control cases in one leaf node
####
		ctfit <- .Call("causalTree",PACKAGE="htetree",
					   ncat = as.integer(cats * !isord),
					   split_Rule = as.integer(split.Rule.int), # tot, ct, fit, tstats, totD, ctD, fitD, tstatsD
					   bucketNum = as.integer(bucketNum), # if == 0, no discrete; else do discrete
					   bucketMax = as.integer(bucketMax),
					   method = as.integer(method.int),  # "anova" not changed yet
					   crossmeth = as.integer(cv.option.num), # tot, ct, fit, tstats
					   crossHonest = as.integer(cv.Honest.num),
					   as.double(unlist(controls)), # control list in rpart
					   minsize, # minsize = min_node_size
					   as.double(propensity),
					   as.integer(xval),
					   as.integer(xgroups),
					   as.double(t(init$y)),
					   X, # X features for model data
					   wt, # for model data
					   treatment, # for model data
					   as.integer(init$numy),
					   as.double(cost),
					   as.double(xvar), # for model daa
					   as.double(split.alpha),
					   as.double(cv.alpha),
					   as.integer(HonestSampleSize),
					   as.double(cv.gamma)
					   )

		nsplit <- nrow(ctfit$isplit) # total number of splits, primary and surrogate
		## total number of categorical splits
		ncat <- if (!is.null(ctfit$csplit)) nrow(ctfit$csplit) else 0L

		if (nsplit == 0L) xval <- 0L # No xvals were done if no splits were found

		numcp <- ncol(ctfit$cptable)
		temp <- if (nrow(ctfit$cptable) == 3L) c("CP", "nsplit", "rel error")
			else c("CP", "nsplit", "rel error", "xerror", "xstd")
		#print (ctfit$cptable)
		dimnames(ctfit$cptable) <- list(temp, 1L:numcp)

		tname <- c("<leaf>", colnames(X))
		# print (tname)
		splits <- matrix(c(ctfit$isplit[, 2:3], ctfit$dsplit), ncol = 5L,
						 dimnames = list(tname[ctfit$isplit[, 1L] + 1L],
										 c("count", "ncat", "improve", "index", "adj")))
		index <- ctfit$inode[, 2L]  # points to the first split for each node

		## Now, make ordered factors look like factors again (a printout choice)
		nadd <- sum(isord[ctfit$isplit[, 1L]])
		if (nadd > 0L) { # number of splits at an ordered factor.
			newc <- matrix(0L, nadd, max(cats))
			cvar <- ctfit$isplit[, 1L]
			indx <- isord[cvar]             # vector of TRUE/FALSE
			cdir <- splits[indx, 2L]        # which direction splits went
			ccut <- floor(splits[indx, 4L]) # cut point
			splits[indx, 2L] <- cats[cvar[indx]] # Now, # of categories instead
			splits[indx, 4L] <- ncat + 1L:nadd # rows to contain the splits

			## Next 4 lines can be done without a loop, but become indecipherable
			for (i in 1L:nadd) {
				newc[i, 1L:(cats[(cvar[indx])[i]])] <- -as.integer(cdir[i])
				newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
			}
			catmat <- if (ncat == 0L) newc
				else {
					## newc may have more cols than existing categorical splits
					## the documentation says that levels which do no exist are '2'
					## and we add 2 later, so record as 0 here.
					cs <- ctfit$csplit
					ncs <- ncol(cs); ncc <- ncol(newc)
					if (ncs < ncc) cs <- cbind(cs, matrix(0L, nrow(cs), ncc - ncs))
					rbind(cs, newc)
				}
			ncat <- ncat + nadd
		} else catmat <- ctfit$csplit

		if (nsplit == 0L) {
			frame <- data.frame(row.names = 1L,
								var = "<leaf>",
								n = ctfit$inode[, 5L],
								wt = ctfit$dnode[, 3L],
								dev = ctfit$dnode[, 1L],
								yval = ctfit$dnode[, 4L],
								complexity = ctfit$dnode[, 2L],
								ncompete = 0L,
								nsurrogate = 0L)
		} else {
			temp <- ifelse(index == 0L, 1L, index)
			svar <- ifelse(index == 0L, 0L, ctfit$isplit[temp, 1L]) # var number
			frame <- data.frame(row.names = ctfit$inode[, 1L],
								var = tname[svar + 1L],
								n = ctfit$inode[, 5L],
								wt = ctfit$dnode[, 3L],
								dev = ctfit$dnode[, 1L],
								yval = ctfit$dnode[, 4L],
								complexity = ctfit$dnode[, 2L],
								ncompete = pmax(0L, ctfit$inode[, 3L] - 1L),
								nsurrogate = ctfit$inode[, 4L])
		}

		if (is.null(init$summary))
			stop("Initialization routine is missing the 'summary' function")
		functions <- if (is.null(init$print)) list(summary = init$summary)
			else list(summary = init$summary, print = init$print)
		if (!is.null(init$text)) functions <- c(functions, list(text = init$text))
		if (method == "user") functions <- c(functions, mlist)

		where <- ctfit$which
		names(where) <- row.names(X)

		ans <- list(frame = frame,
					where = where,
					call = Call, terms = Terms,
					cptable = t(ctfit$cptable),
					method = method,
					control = controls,
					functions = functions,
					numresp = init$numresp)
		if (nsplit) ans$splits = splits
		if (ncat > 0L) ans$csplit <- catmat + 2L
		if (nsplit) ans$variable.importance <- importance(ans)

		if (y) ans$y <- Y
		if (x) {
			ans$x <- X
			ans$wt <- wt
		}
		ans$ordered <- isord
		if (!is.null(attr(m, "na.action"))) ans$na.action <- attr(m, "na.action")
		if (!is.null(xlevels)) attr(ans, "xlevels") <- xlevels
		if (method == "class") attr(ans, "ylevels") <- init$ylevels
		class(ans) <- "rpart"
		if(ncol(ans$cptable) >= 4) {
			ans$cptable[,4]  <- ans$cptable[,4] / ans$cptable[1, 4]
		}
		ans
}
