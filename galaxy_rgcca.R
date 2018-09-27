# Author: Etienne CAMENEN
# Date: 2018
# Institute: ICM - Institut du Cerveau et de la Moelle ?pini?re (Paris, FRANCE),
# Institut Fran?ais de Bioinformatique (IFB), Centre national de la recherche scientifique (CNRS)
# Contact: iconics@icm-institute.org
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
# EDAM topic: omics, medecine, mathematics
#
# Abstract: A user-friendly multi-blocks analysis (Regularized Generalized Canonical Correlation Analysis, RGCCA)
# with all default settings predefined. Produce two figures to help clinicians to identify biomarkers: 
# samples and variables projected on the two first component of the multi-block analysis.

rm(list=ls())

################################
#          File
################################

getFileName = function(fi)  {
  #get prefix part from a file
  fo = unlist(strsplit(fi, '/'))
  fo = fo[length(fo)]
  fo = unlist(strsplit(fo, '[.]'))[1]
}

loadData = function(fi, fo=fi, row.names=NULL, h=F){
  #create a dataset object from a file loading
  #fi: input file name
  #fo: dataset object name
  data = as.matrix(read.table(fi, sep = SEPARATOR, h = h, row.names = row.names))
  assign(fo, data, .GlobalEnv)
  #TODO: catch warning missing \n at the end of the file
}

save = function(f, p){
  pdf(f, width=10, height=8)
  plot(p)
  suprLog = dev.off()
}

setBlocks = function(){
  #create a list object of blocks from files loading
  
  #remove white space
  opt$datasets = gsub(" ", "", opt$datasets)
  #split by ,
  blocksName = unlist(strsplit(opt$datasets, ","))
  
  #load each dataset
  blocks = list()
  for (i in 1:length(blocksName)){
    fi = blocksName[i]
    checkFile(fi)
    fo = getFileName(fi)
    loadData(fi, fo, 1, T)
    blocks[[fo]] = get(fo)
    if (NCOL(blocks[[fo]]) ==0) stop(paste(fo, "block file has an only-column. Check the --separator [by default: 1 for tabulation].\n"), call.=FALSE)
  }
  blocks[["Superblock"]] = Reduce(cbind, blocks)
  
  return(blocks)
}


setConnection = function(){
  #default settings of connection_matrix matrix
  if(is.null(opt$connection)){
    seq = 1:(length(blocks)-1)
    connection_matrix = matrix(0,length(blocks),length(blocks))
    connection_matrix[length(blocks), seq] <- 1 -> connection_matrix[seq, length(blocks)]
  }else{
    loadData(opt$connection, "connection_matrix", h=F)
  }
  return(connection_matrix)
}

setResponse = function(){
  #create a dataset object from a file loading containg the response
  
  if("response" %in% names(opt)){
    loadData(opt$response, "response", 1, F)
    #TODO: check n1  = n2 = ...
    if(isTRUE(DISJONCTIF)) response = factor(apply("Response", 1, which.max))
    return (response)
  }else{
    return ( rep("black", NROW(blocks[[1]])) )
  }
}

################################
#          Graphic
################################

circleFun <- function(center = c(0,0), diameter = 2, npoints = 100){
  #creates x,y coordinates for a circle
  
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

printAxis = function (n)
  #n: number of the axis
  paste("Axis ", n, " (", round(rgcca$AVE$AVE_X[[length(blocks)]][n] * 100 , 1),"%)", sep="")

theme_perso = function() {
  theme(
    legend.text = element_text(size = 13),
    legend.title = element_text(face="bold.italic", size=16),
    plot.title = element_text(size = 25, face = "bold", hjust=0.5, margin = margin(0,0,20,0))
  )
}

plotSpace = function (df, title, response, name_group, comp1, comp2){
  #plot settings for projection of points in a bi-dimensional space
  ggplot(df, aes(df[,1], df[,2], colour = response)) + 
  theme_classic() +
  geom_vline(xintercept = 0, col="grey", linetype="dashed", size=1) + 
  geom_hline(yintercept = 0, col="grey", linetype="dashed", size=1) + 
  labs ( title = paste(title, "space"),
         x = printAxis(comp1), 
         y = printAxis(comp2),
         color = name_group) +
  geom_text_repel(aes(colour = response, label= rownames(df)), size = 3, force=2) +
  scale_y_continuous(breaks=NULL) +
  scale_x_continuous(breaks=NULL) +
  theme_perso() +
  theme(
    axis.text = element_blank(),
    axis.title.y = element_text(face=AXIS_FONT, margin = margin(0,20,0,0), size=AXIS_TITLE_SIZE),
    axis.title.x = element_text(face=AXIS_FONT, margin = margin(20,0,0,0), size=AXIS_TITLE_SIZE)
    )
  #+ stat_ellipse()
  #TODO: if NB_VAR > X
}

#TODO: convert coef into [-1,1]
plot_biomarkers = function(df, comp, n){
  
  df = data.frame(df[order(abs(df[,comp]), decreasing = TRUE),], order = nrow(df):1)
  color2=df$color; levels(color2)=hue_pal()(length(blocks)-1)
  if(NROW(df) >= n) df = df[1:n,]
  
  ggplot(df, aes(order, df[,comp], fill = color)) +
  geom_hline(yintercept = c(-.5,.5), col="grey", linetype="dotted", size=1) + 
  geom_hline(yintercept = 0, col="grey", size=1) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_x_continuous(breaks=df$order, labels=rownames(df)) +
  scale_y_continuous(breaks=seq(-1,1,.5), limits = c(-1,1)) +
  labs(
    title= "Variable weights", 
    subtitle=printAxis(comp),
    x = "", y = "",
    fill = "Blocks") +
  theme_classic() +
  theme_perso() +
  theme(
        axis.text.y = element_text(size = AXIS_TEXT_SIZE, face=AXIS_FONT, color=as.character(color2)),
        axis.text.x = element_text(size = AXIS_TEXT_SIZE, face=AXIS_FONT, color="darkgrey"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, size = 16, face="italic"))
  
}

################################
#          Arguments
################################

#TODO: remove default files
getArgs = function(){
  option_list = list(
    make_option(c("-d", "--datasets"), type="character", metavar="character", default = "data/agriculture.tsv,data/industry.tsv,data/politic.tsv", help="Bloc files name"),
    make_option(c("-c", "--connection"), type="character", metavar="character", help="Connection file name"),
    make_option(c("-r", "--response"), type="character", metavar="character", help="Response file name"),
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

checkFile = function (o){
  # o: one argument from the list of arguments
  if(!file.exists(o)){
    stop(paste(o, " file does not exist\n", sep=""), call.=FALSE)
  }
}

#Check the validity of the arguments 
#Inputs:
# a: arguments (optionParser object)
checkArg = function(a){
  opt = parse_args(a)
  
  print(opt$datasets)
  if(is.null(opt$datasets)) stop(paste("--datasets is required\n", sep=""), call.=FALSE)
  
  if (is.null(opt$scheme)) opt$scheme = "factorial"
  else if ((opt$scheme < 1) || (opt$scheme > 4)){
    stop("--scheme must be comprise between 1 and 4 [by default: 3].\n", call.=FALSE)
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

################################
#            MAIN
################################



#Loading librairies
librairies = c("RGCCA", "ggplot2", "ggrepel", "optparse", "scales")
for (l in librairies){
  if (! (l %in% installed.packages()[,"Package"])) install.packages(l, repos = "http://cran.us.r-project.org", quiet = T)
  library(l, character.only = TRUE)
}

#Global settings
SCALE = T
VERBOSE = F
TAU = "optimal"
DISJONCTIF = F
COMP1 = 1
COMP2 = 2
AXIS_TITLE_SIZE = 19
AXIS_TEXT_SIZE = 10
AXIS_FONT = "italic"

#Get arguments
args = getArgs()
tryCatch({
  opt = checkArg(args)
}, error = function(e) {
  #print_help(args)
  stop(e[[1]], call.=FALSE)
})



blocks = setBlocks()
print(names(blocks))
connection_matrix = setConnection()
response = setResponse()
NB_COMP = sapply(blocks, NCOL)
# TODO: Error in rgcca(blocks, connection_matrix, tau = TAU, scheme = scheme, ncomp = rep(NB_COMP,  : 
#                                                                     For each block, choose a number of components smaller than the number of variables!

rgcca = rgcca(blocks,
              connection_matrix,
              tau = TAU,
              scheme = opt$scheme,
              ncomp = NB_COMP,
              scale = SCALE,
              verbose = VERBOSE)
#ncomp = rep(NB_COMP, length(blocks))
#TODO: catch Error in connection_matrix * h(cov2(Y, bias = bias)) : non-conformable arrays
#message: Number of row/column of connection matrix doesn't match with the number of blocks.


# Samples common space
samples = data.frame(rgcca$Y[[length(blocks)]])
samplesSpace = plotSpace(samples, "Samples", response, "Response", COMP1, COMP2)
save(opt$output1, samplesSpace)

#attribution of block ID to each corresponding variable
blocks_variables = rep( names(blocks)[-length(blocks)] , sapply(blocks[1:(length(blocks)-1)], NCOL))

# Variables common space
variables =  data.frame( 
 #correlation matrix with superblock for each variables and each component selected
 sapply ( c(COMP1:COMP2), function(x) cor( blocks[["Superblock"]], rgcca$Y[[length(blocks)]][, x] ) ) , 
 blocks_variables,
 row.names = colnames(blocks[["Superblock"]])
)

variablesSpace = plotSpace(variables, "Variables", variables[,3], "Blocks", COMP1, COMP2) + 
  geom_path(aes(x,y), data=circleFun(), col="grey", size=1) + 
  geom_path(aes(x,y), data=circleFun()/2, col="grey", size=1, lty=2)
save(opt$output2, variablesSpace)

# Biomarkers plot
biomarkers = data.frame(rgcca$a[[length(blocks)]], color=blocks_variables)
best_biomarkers = plot_biomarkers(biomarkers, 1, 10)
save(opt$output3, best_biomarkers)