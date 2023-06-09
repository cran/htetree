% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveFiles.R
\name{saveFiles}
\alias{saveFiles}
\title{Save Necessary Files to Run Shiny App}
\usage{
saveFiles(
  model,
  data,
  outcomevariable,
  treatment_indicator,
  propensity_score = "",
  filePath = ""
)
}
\arguments{
\item{model}{a tree model constructed by \code{hte_causalTree, hte_matchinleaves,
or hte_ipw}.}

\item{data}{a data frame containing the variables in the model.}

\item{outcomevariable}{a character representing the column name
of the outcome variable.}

\item{treatment_indicator}{a character representing the column name
for the treatment variable in the causal setup.}

\item{propensity_score}{a character representing the column name of
the propensity score.}

\item{filePath}{a character string representing the path name
to save the files temporarily.}
}
\value{
No return value. It is used to save necessary files temporarily
to run Shiny App.
}
\description{
This function is to save files necessary to run Shiny app to
visualize causal tree and the estimated heterogeneous treatment
effects in an interactive way.
}
