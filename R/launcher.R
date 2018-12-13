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
# with all default settings predefined. Produce two figures to help clinicians to identify biomarkers:
# samples and variables projected on the two first component of the multi-block analysis.

rm(list=ls())

##################
#     Arguments
##################

getArgs = function(){
  option_list = list(
    make_option(c("-d", "--datasets"), type="character", metavar="character", help="Path of the blocks", default = "data2/Clinique.tsv,data2/Lipidomique.tsv,data2/Transcriptomique.tsv,data2/Imagerie.tsv,data2/Metabolomique.tsv"),
    make_option(c("-w", "--directory"), type="character", metavar="character", help="Path of the scripts directory (for Galaxy)", default="."),
    make_option(c("-c", "--connection"), type="character", metavar="character", help="Connection file path"),
    make_option(c("-r", "--response"), type="character", metavar="character", help="Response file path"),
    make_option(c("-n", "--names"), type="character", metavar="character", help="Names of the blocks [default: filename]"),
    make_option(c("-H", "--header"), type="logical", action="store_false", help="DO NOT consider first row as header of columns"),
    make_option(c("-s", "--separator"), type="integer", metavar="integer", default=1,
                help="Type of separator [default: tabulation] (1: Tabulation, 2: Semicolon, 3: Comma"),
    make_option(c("-g", "--scheme"), type="integer", metavar="integer", default=2,
                help="Scheme function g(x) [default: x^2] (1: x, 2: x^2, 3: |x|, 4: x^4"),
    make_option(c( "--output1"), type="character", metavar="character", default="samples_space.pdf",
                help="Variables space file name [default: %default]"),
    make_option(c( "--output2"), type="character", metavar="character", default="variables_space.pdf",
                help="Sample space file name [default: %default]"),
    make_option(c( "--output3"), type="character", metavar="character", default="best_biomarkers.pdf",
                help="Best biomarkers file name [default: %default]")
  )
  args = commandArgs(trailingOnly=T)
  return (OptionParser(option_list=option_list))
}



#Check the validity of the arguments
#Inputs:
# a: arguments (optionParser object)
checkArg = function(a){
  opt = parse_args(a)

  if(is.null(opt$datasets)) stop(paste("--datasets is required\n", sep=""), call.=FALSE)

  if (is.null(opt$scheme)) opt$scheme = "factorial"
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
    assign("SEPARATOR", separators[opt$separator], .GlobalEnv)
  }

  FILES = c("connection", "response")
  for (o in FILES)
    if(!is.null(opt[[o]])) checkFile(opt[[o]])

  return (opt)
}

checkFile = function (f){
  # o: one argument from the list of arguments
  if(!file.exists(f)){
    stop(paste(f, " file does not exist\n", sep=""), call.=FALSE)
  }
}

##################
#     Main
##################

# Pre-requisite: for xlsx inputs, java must be installed
# Under linux: sudo apt-get install default-jre default-jdk && sudo R CMD javareconf

#Loading librairies
librairies = c("RGCCA", "ggplot2", "ggrepel", "optparse", "scales", "cluster", "xlsx")
for (l in librairies) {
  if (!(l %in% installed.packages()[, "Package"]))
    install.packages(l, repos = "http://cran.us.r-project.org", quiet = T)
  library(l, character.only = TRUE)
}

#Get arguments
args = getArgs()
tryCatch({
  opt = checkArg(args)
}, error = function(e) {
  #print_help(args)
  stop(e[[1]], call.=FALSE)
})

#Global settings
SCALE = T
VERBOSE = F
TAU = "optimal"
DISJUNCTIVE = F
COMP1 = 1
COMP2 = 2
AXIS_TITLE_SIZE = 19
AXIS_TEXT_SIZE = 8
PCH_TEXT_SIZE = 2
AXIS_FONT = "italic"
MAX_CLUSTERS = 10
COLOR_SAMPLES_DEF = "#000099"
HEADER = !("header" %in% names(opt))
MSG_HEADER = " Possible mistake: header parameter is disabled, check if the file does'nt have one."
NB_MARK = 100
SUPERBLOCK = T

setwd(opt$directory)
source("R/parsing.R")
source("R/plot.R")

blocks = setBlocks()
connection_matrix = setConnection()
response = setResponse()
NB_COMP = 2
ncomp = rep(NB_COMP, length(blocks))
#sapply(blocks, NCOL)
# TODO: Error in rgcca(blocks, connection_matrix, tau = TAU, scheme = scheme, ncomp = rep(NB_COMP,  :
#                                                                     For each block, choose a number of components smaller than the number of variables!

getColumnSameVal = function(list_m)
  lapply(1:length(list_m), function (x) which( apply(list_m[[x]], 2, sd ) == 0 ))
#getColumnSameVal(blocks)

rsgcca.res = sgcca(A = blocks,
              C = connection_matrix,
              scheme = opt$scheme,
              ncomp = ncomp,
              scale = SCALE,
              verbose = VERBOSE)

# Samples common space
( samplesSpace = plotSamplesSpace(rsgcca.res, COMP1, COMP2) )
plotSamplesSpace(rsgcca.res, COMP1, COMP2, 1)
save(opt$output1, samplesSpace)

# Variables common space
( variablesSpace = plotVariablesSpace(rsgcca.res, COMP1, COMP2) )
plotVariablesSpace(rsgcca.res, COMP1, COMP2, 2)
save(opt$output2, variablesSpace)

# Biomarkers plot
( best_biomarkers = plot_biomarkers(rsgcca.res, COMP1, NB_MARK) )
plot_biomarkers(rsgcca.res, COMP1, NB_MARK, 2)
save(opt$output3, best_biomarkers)

plotAVE(rsgcca.res, COMP1)
