# Author: Etienne CAMENEN
# Date: 2018
# Institute: ICM - Institut du Cerveau et de la Moelle epiniere (Paris, FRANCE),
# Institut Francais de Bioinformatique (IFB), Centre national de la recherche scientifique (CNRS)
# Contact: iconics@icm-institute.org
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
# EDAM topic: omics, medecine, mathematics
#
# Abstract: A user-friendly multi-blocks analysis (Regularized Generalized Canonical Correlation Analysis, RGCCA)
# with all default settings predefined. Produce four figures to help clinicians to identify fingerprint:
# the samples and the variables projected on the two first component of the multi-block analysis, the histograms
# of the most explicative variables and the explained variance for each blocks.


rm(list=ls())

##################
#     Arguments
##################

# Parse the arguments from a command line launch
getArgs = function(){
  option_list = list(
    make_option(c("-d", "--datasets"), type="character", metavar="character", help="Path of the block files", default = opt[15]),
    make_option(c("-w", "--directory"), type="character", metavar="character", help="Path of the scripts directory (for Galaxy)", default=opt[1]),
    make_option(c("-c", "--connection"), type="character", metavar="character", help="Path of the connection file"),
    make_option(c("-r", "--response"), type="character", metavar="character", help="Path of the response file"),
    make_option(c("--names"), type="character", metavar="character", help="Names of the blocks [default: filename]"),
    make_option(c("-H", "--header"), type="logical", action="store_false", help="DO NOT consider first row as header of columns"),
    make_option(c("--separator"), type="integer", metavar="integer", default=1,
                help="Type of separator (1: Tabulation, 2: Semicolon, 3: Comma) [default: tabulation]"),
    make_option(c("-t", "--tau"), type="character", metavar="character/double", default=opt[4],
                help="Tau parameter for RGCCA (eihter a double between 0 and 1 or the character 'optimal' for an automatic setting)"),
    make_option(c("-g", "--scheme"), type="integer", metavar="integer", default=2,
                help="Scheme function g(x) (1: x, 2: x^2, 3: |x|, 4: x^4) [default: x^2]"),
    make_option(c("--scale"),  type="logical", action="store_false",
                help="DO NOT scale the blocks"),
    make_option(c("--superblock"),  type="logical", action="store_false",
                help="DO NOT use a superblock"),
    make_option(c("--init"),  type="integer", metavar="integer", default=1,
                help="Initialization mode (1: Singular Value Decompostion , 2: random) [default: SVD]"),
    make_option(c("--bias"),  type="logical", action="store_false",
                help="Unbiased estimator of the var/conv"),
    make_option(c("--ncomp"),  type="integer", metavar="integer", default=opt[6],
                help="Number of components in the analysis"),
    make_option(c("--block"),  type="integer", metavar="integer", default=opt[7],
                help="Number of the block shown in the graphics (0: the superblock or, if not, the last, 1: the fist one, 2: the 2nd, etc.) [default: the last one]"),
    make_option(c("--compx"),  type="integer", metavar="integer", default=opt[8],
                help="X-axis for biplots and component for histograms"),
    make_option(c("--compy"),  type="integer", metavar="integer", default=opt[9],
                help="Y-axis for biplots"),
    make_option(c("--nmark"),  type="integer", metavar="integer", default=opt[10],
                help="Number maximum of bioimarkers in the fingerprint"),
    make_option(c( "--output1"), type="character", metavar="character", default=opt[11],
                help="Variables space file name [default: %default]"),
    make_option(c( "--output2"), type="character", metavar="character", default=opt[12],
                help="Sample space file name [default: %default]"),
    make_option(c( "--output3"), type="character", metavar="character", default=opt[13],
                help="Best fingerprint file name [default: %default]"),
    make_option(c( "--output4"), type="character", metavar="character", default=opt[14],
                help="AVE plot file name [default: %default]")
  )
  args = commandArgs(trailingOnly=T)
  return (OptionParser(option_list=option_list))
}

checkFile = function (f){
  # Check the existence of a path
  # f: A character giving the path of a file

  if(!file.exists(f)){
    stop(paste(f, " file does not exist\n", sep=""), call.=FALSE)
  }
}

# Check the validity of the arguments
# opt : an optionParser object
checkArg = function(opt){

  if(is.null(opt$datasets))
    stop(paste("--datasets is required\n", sep=""), call.=FALSE)

  if (is.null(opt$scheme))
    opt$scheme = "factorial"
  else if ((opt$scheme < 1) || (opt$scheme > 4)){
    stop("--scheme must be comprise between 1 and 4 [by default: 2].\n", call.=FALSE)
  }else{
    schemes = c("horst", "factorial", "centroid")
    if (opt$scheme == 4)
      opt$scheme = function(x) x^4
    else
      opt$scheme = schemes[opt$scheme]
  }

  if ((opt$separator < 1) || (opt$separator > 3)){
    stop("--separator must be comprise between 1 and 2 (1: Tabulation, 2: Semicolon, 3: Comma) [by default: 2].\n", call.=FALSE)
  }else{
    separators = c('\t', ';', ',')
    opt$separator = separators[opt$separator]
  }

  if ((opt$init < 1) || (opt$init > 2)){
    stop("--init must be comprise between 1 and 2 (1: Singular Value Decompostion , 2: random) [by default: SVD].\n", call.=FALSE)
  }else{
    opt$init = ifelse(opt$init == 1, "svd", "random")
  }

  FILES = c("connection", "response")
  for (o in FILES)
    if(!is.null(opt[[o]]))
      checkFile(opt[[o]])

  return (opt)
}

# Check the validity of the arguments after loading the blocks
# opt : an optionParser object
# blocks : a list of matrix
postCheckArg = function(opt, blocks){
  if ((opt$ncomp < 2) || (opt$ncomp > min(sapply(blocks, NCOL)))){
    stop("--ncomp must be comprise between 2 and ", min(sapply(blocks, NCOL)) ," (the minimum number of variables among the whole blocks).\n", call.=FALSE)
  }

  checkComp = function (x){
    if ((opt[[x]] < 1) || (opt[[x]] > opt$ncomp )){
      stop(paste("--", x, " must be comprise between 2 and ", opt$ncomp ," (the number of component selected).\n ", sep=""), call.=FALSE)
    }
  }
  out = sapply(c("compx", "compy"), function(x) checkComp (x))

  MSG = "--tau must be comprise between 0 and 1 or must corresponds to the character 'optimal' for automatic setting.\n"
  if (opt$tau != "optimal"){
    tryCatch({
      opt$tau = as.double(gsub(",", ".",  opt$tau))
      if((opt$tau < 0) || (opt$tau > 1))
        stop(MSG, call.=FALSE)
      else
        opt$tau = rep(opt$tau, length(blocks))
    }, warning = function(w) {
      stop(MSG, call.=FALSE)
    })
  }

  if(opt$block > length(blocks))
    stop(paste("--block must be lower than ", length(blocks), " (the maximum number of blocks).\n", sep=""), call.=FALSE)
  else if(opt$block == 0)
    opt$block = length(blocks)

  # if (opt$nmark > NCOL(blocks[[opt$block]])){
  #   stop(paste("--nmark must be lower than ", NCOL(blocks[[opt$block]]) ," (the maximum number of columns among the block selected).\n", sep=""), call.=FALSE)
  # }

  return (opt)
}

#' Launch a Shiny application for S/RGCCA
#' @export
runShiny <- function()
  shiny::runApp("inst/shiny")

##################
#     Main
##################

# Pre-requisite: for xlsx inputs, java must be installed
# Under linux: sudo apt-get install default-jre default-jdk && sudo R CMD javareconf

#Loading librairies
librairies = c("RGCCA", "ggplot2", "optparse", "scales", "xlsx")
for (l in librairies) {
  if (!(l %in% installed.packages()[, "Package"]))
    install.packages(l, repos = "http://cran.us.r-project.org", quiet = T)
  library(l, character.only = TRUE)
}

# Get arguments : R packaging install, need an opt variable with associated arguments
opt = list(directory = ".",
           separator = "\t",
           scheme = "factorial",
           tau = "optimal",
           init = "svd",
           ncomp = 2,
           block = 0,
           compx = 1,
           compy = 2,
           nmark = 100,
           output1 = "samples_plot.pdf",
           output2 = "corcircle.pdf",
           output3 = "fingerprint.pdf",
           output4 = "ave.pdf",
           datasets="data2/Clinique.tsv,data2/Lipidomique.tsv,data2/Transcriptomique.tsv,data2/Imagerie.tsv,data2/Metabolomique.tsv")

tryCatch({
  opt = parse_args(getArgs())
  opt = checkArg(opt)
}, error = function(e) {
  if (length(grep("nextArg", e[[1]])) != 1)
    stop(e[[1]], call.=FALSE)
})

setwd(opt$directory)
source("R/parsing.R")
source("R/plot.R")

# Global settings
opt$header = !("header" %in% names(opt))
opt$superblock = !("superblock" %in% names(opt))
opt$bias = !("bias" %in% names(opt))
opt$scale = !("scale" %in% names(opt))
VERBOSE = F

blocks = setBlocks(opt$superblock, opt$datasets, opt$names, opt$separator, opt$header)
opt = postCheckArg(opt, blocks)
connection = setConnection(blocks, opt$connection, opt$separator)
response = setResponse(blocks, opt$response, opt$separator, opt$header)
# ncomp = sapply(blocks, NCOL)
# TODO: Error in rgcca(blocks, connection_matrix, tau = TAU, scheme = scheme, ncomp = rep(NB_COMP,  :
#                                                                     For each block, choose a number of components smaller than the number of variables!

getColumnSdNull = function(list_m)
  lapply(list_m, function (x) which( apply(x, 2, sd ) == 0 ))
# TODO: getColumnSdNull(blocks)

rgcca.res = rgcca(A = blocks,
              C = connection,
              scheme = opt$scheme,
              ncomp = rep(opt$ncomp, length(blocks)),
              scale = opt$scale,
              tau = opt$tau,
              verbose = VERBOSE,
              init = opt$init,
              bias = opt$bias)

names(rgcca.res$a) = names(blocks)

# Samples common space
( samples_plot = plotSamplesSpace(rgcca.res, response, opt$compx, opt$compy, opt$block) )
plotSamplesSpace(rgcca.res, response, opt$compx, opt$compy, 1)
savePlot(opt$output1, samples_plot)

# Variables common space
( corcircle = plotVariablesSpace(rgcca.res, blocks, opt$compx, opt$compy, opt$superblock, opt$block) )
plotVariablesSpace(rgcca.res, blocks, opt$compx, opt$compy, opt$superblock, 1)
savePlot(opt$output2, corcircle)

# Fingerprint plot
( fingerprint = plotFingerprint(rgcca.res, opt$compx, opt$superblock, opt$nmark, opt$block) )
plotFingerprint(rgcca.res, opt$compx, opt$superblock, opt$nmark, 2)
savePlot(opt$output3, fingerprint)

# Average Variance Explained
ave = plotAVE(rgcca.res, opt$compx)
savePlot(opt$output4, ave)
