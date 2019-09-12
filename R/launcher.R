# Author: Etienne CAMENEN
# Date: 2019
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

#' @import RGCCA
#' @import ggplot2
#' @importFrom grDevices dev.off rgb colorRamp pdf colorRampPalette
#' @importFrom graphics plot
#' @importFrom stats cor quantile runif sd na.omit p.adjust pnorm qnorm weights
#' @importFrom utils read.table write.table packageVersion installed.packages head
#' @importFrom scales hue_pal
#' @importFrom optparse OptionParser make_option parse_args
#' @importFrom plotly layout ggplotly style plotly_build %>% plot_ly add_trace
#' @importFrom visNetwork visNetwork visNodes visEdges
#' @importFrom igraph graph_from_data_frame V<- E<-
#' @importFrom methods is

rm(list = ls())
graphics.off()

########## Arguments ##########

# Parse the arguments from a command line launch
getArgs <- function() {
    option_list <- list(
        # File parameters
        make_option(
            opt_str = c("-d", "--datasets"),
            type = "character",
            metavar = "path list",
            default = opt[20],
            help = "List of comma-separated file paths corresponding to the
            blocks to be analyzed (one per block and without spaces between
            them; e.g., path/file1.txt,path/file2.txt) [required]"
        ),
        make_option(
            opt_str = c("-w", "--directory"),
            type = "character",
            metavar = "path",
            default = opt[1],
            help = "Path of the root folder containing the R/ (e.g. for Galaxy)
            [default: the current one]"
        ),
        make_option(
            opt_str = c("-c", "--connection"),
            type = "character",
            metavar = "path",
            help = "Path of the file defining the connections between the blocks
            [if not used, activates the superblock mode]"
        ),
        make_option(
            opt_str = "--group",
            type = "character",
            metavar = "path",
            help = "Path of the file coloring the individuals in the ad hoc
            plot"
        ),
        make_option(
            opt_str = c("-r", "--response"),
            type = "integer",
            metavar = "integer",
            help = "Position of the response file for the supervised mode within
            the block path list [actives the supervised mode]"
        ),
        make_option(
            opt_str = "--names",
            type = "character",
            metavar = "character list",
            help = "List of comma-separated block names to rename them (one per
            block; without spaces between them) [default: the block file names]"
        ),
        make_option(
            opt_str = c("-H", "--header"),
            type = "logical",
            action = "store_false",
            help = "DO NOT consider the first row as the column header"
        ),
        make_option(
            opt_str = "--separator",
            type = "integer",
            metavar = "integer",
            default = 1,
            help = "Character used to separate columns (1: tabulation,
            2: semicolon, 3: comma) [default: %default]"
        ),
        # Analysis parameter
        make_option(
            opt_str = "--type",
            type = "character",
            metavar = "character",
            default = opt[3],
            help = "Type of analysis [default: %default] (among: rgcca, pca,
            cca, gcca, cpca-w, hpca, maxbet-b, maxbet, maxdiff-b, maxdiff,
            maxvar-a, maxvar-b, maxvar, niles, r-maxvar, rcon-pca, ridge-gca,
            sabscor, ssqcor, ssqcor, ssqcov-1, ssqcov-2, ssqcov, sum-pca,
            sumcor, sumcov-1, sumcov-2, sumcov)"
        ),
        make_option(
            opt_str = "--ncomp",
            type = "character",
            metavar = "integer list",
            default = opt[4],
            help = "Number of components in the analysis for each block
            [default: %default]. The number should be greater than 1 and lower
            than the minimum number of variables among the blocks. It can be a
            single values or a comma-separated list (e.g 2,2,3,2)."
        ),
        make_option(
            opt_str = "--tau",
            type = "character",
            metavar = "float list",
            default = opt[5],
            help = "A regularization parameter for each block (i.e., tau)
            [default: %default]. Tau varies from 0 (maximizing the correlation)
            to 1 (maximizing the covariance). For SGCCA, tau is automatically
            set to 1. A shrinkage parameter can be defined instead for
            automatic variable selection, varying from the square root of the
            variable number (the fewest selected variables) to 1 (all the
            variables are included). It can be a single value or a
            comma-separated list (e.g. 0,1,0.75,1)."
        ),
        make_option(
            opt_str = "--scheme",
            type = "integer",
            metavar = "integer",
            default = 2,
            help = "Link (i.e. scheme) function for covariance maximization
            (1: x, 2: x^2, 3: |x|, 4: x^4) [default: %default]. Only, the x
            function penalizes structural negative correlation. The x^4
            function discriminates more strongly the blocks than the x^2 one."
        ),
        make_option(
            opt_str = "--scale",
            type = "logical",
            action = "store_false",
            help = "DO NOT scale the blocks (i.e., a data centering step is
            always performed). Otherwise, each block is normalised and divided
            by the squareroot of its number of variables."
        ),
        make_option(
            opt_str = "--superblock",
            type = "logical",
            action = "store_false",
            help = "DO NOT use a superblock (i.e. a concatenation of all the
            blocks to visualize them all together in a consensus space). In
            this case, all blocks are assumed to be connected or a connection
            file could be used."
        ),
        # Graphical parameters
        make_option(
            opt_str = "--text",
            type = "logical",
            action = "store_false",
            help = "Display the name of the points instead of shapes when
            plotting"
        ),
        make_option(
            opt_str = "--block",
            type = "integer",
            metavar = "integer",
            default = opt[8],
            help = "Position in the path list of the plotted block (0: the
            superblock or, if not activated, the last one, 1: the fist one,
            2: the 2nd, etc.)[default: the last one]"
        ),
        make_option(
            opt_str = "--block_y",
            type = "integer",
            metavar = "integer",
            help = "Position in the path list of the plotted block for the
            Y-axis in the individual plot (0: the superblock or, if not
            activated, the last one, 1: the fist one, 2: the 2nd, etc.)
            [default: the last one]"
        ),
        make_option(
            opt_str = "--compx",
            type = "integer",
            metavar = "integer",
            default = opt[9],
            help = "Component used in the X-axis for biplots and the only
            component used for histograms [default: %default] (should not be
            greater than the number of components of the analysis)"
        ),
        make_option(
            opt_str = "--compy",
            type = "integer",
            metavar = "integer",
            default = opt[10],
            help = "Component used in the Y-axis for biplots
            [default: %default] (should not be greater than the number of
            components of the analysis)"
        ),
        make_option(
            opt_str = "--nmark",
            type = "integer",
            metavar = "integer",
            default = opt[11],
            help = "Number maximum of top variables in ad hoc plot
            [default: %default]"
        ),
        # output parameters
        make_option(
            opt_str = "--o1",
            type = "character",
            metavar = "path",
            default = opt[12],
            help = "Path for the variable plot [default: %default]"
        ),
        make_option(
            opt_str = "--o2",
            type = "character",
            metavar = "path",
            default = opt[13],
            help = "Path for the individual plot [default: %default]"
        ),
        make_option(
            opt_str = "--o3",
            type = "character",
            metavar = "path",
            default = opt[14],
            help = "Path for the top variables plot [default: %default]"
        ),
        make_option(
            opt_str = "--o4",
            type = "character",
            metavar = "path",
            default = opt[15],
            help = "Path for the explained variance plot [default: %default]"
        ),
        make_option(
            opt_str = "--o5",
            type = "character",
            metavar = "path",
            default = opt[16],
            help = "Path for the design plot [default: %default]"
        ),
        make_option(
            opt_str = "--o6",
            type = "character",
            metavar = "path",
            default = opt[17],
            help = "Path for the individual table [default: %default]"
        ),
        make_option(
            opt_str = "--o7",
            type = "character",
            metavar = "path",
            default = opt[18],
            help = "Path for the variable table [default: %default]"
        ),
        make_option(
            opt_str = "--o8",
            type = "character",
            metavar = "path",
            default = opt[19],
            help = "Path for the analysis results in RData [default: %default]"
        )
    )
    return(OptionParser(option_list = option_list))
}

checkFile <- function(f) {
    # Check the existence of a path f: A character giving the path of a file
    
    if (!file.exists(f)) {
        stop(paste0(f, " file does not exist."), exit_code = 120)
    }
}

checkArg <- function(opt) {
        # Check the validity of the arguments opt : an optionParser object
        
        if (is.null(opt$datasets))
            stop(paste0("--datasets is required."), exit_code = 121)
    
    if (is.null(opt$scheme))
        opt$scheme <- "factorial"
    else if (!opt$scheme %in% seq_len(4)) {
        stop(
            paste0(
                "--scheme must be comprise between 1 and 4 [by default: 2], not be equal to ",
                opt$scheme,
                "."
            ),
            exit_code = 122
        )
    } else {
        schemes <- c("horst", "factorial", "centroid")
        if (opt$scheme == 4)
            opt$scheme <- function(x) x ^ 4
        else
            opt$scheme <- schemes[opt$scheme]
    }
    
    if (!opt$separator %in% seq_len(3)) {
        stop(
            paste0(
                "--separator must be comprise between 1 and 3 (1: Tabulation, 2: Semicolon, 3: Comma) [by default: 2], not be equal to ",
                opt$separator,
                "."
            ),
            exit_code = 123
        )
    } else {
        separators <- c("\t", ";", ",")
        opt$separator <- separators[opt$separator]
    }
    
    # if (! opt$init %in% 1:2 )
    # stop(paste0('--init must be 1 or 2 (1: Singular Value
    # Decompostion , 2: random) [by default: 1], not ', opt$init, '.'),
    #  exit_code = 124)
    # else
    # opt$init <- ifelse(opt$init == 1, 'svd', 'random')
    
    
    FILES <- c("connection", "group")
    for (o in FILES)
        if (!is.null(opt[[o]]))
            checkFile(opt[[o]])
    
    checkInteger("nmark")
    if (opt$nmark < 2)
        stop(paste0("--nmark must be upper than 2, not be equal to ",
            opt$nmark,
            "."),
            exit_code = 135)
    
    return(opt)
}

checkArgSize <- function(blocks, x, y) {
    if (length(x) != length(blocks))
        stop(
            paste0(
                "--",
                y,
                " list must have the same size (",
                length(x),
                ") than the number of blocks (",
                length(blocks),
                ")."
            ),
            exit_code = 130
        )
    else
        return(TRUE)
}

postCheckArg <- function(opt, blocks) {
    # Check the validity of the arguments after loading the blocks opt : an
    # optionParser object blocks : a list of matrix
    
    if (!is.null(opt$names))
        checkArgSize(
            blocks,
            strsplit(gsub(" ", "", opt$names), ",")[[1]], "names"
        )
    
    opt <- select_type(blocks, opt)
    
    if (opt$superblock | opt$type == "pca")
        blocks <- c(blocks, list(Reduce(cbind, blocks)))
    
    for (x in c("block", "block_y", "response")) {
        if (!is.null(opt[[x]])) {
            if (opt[[x]] > length(blocks))
                stop(
                    paste0(
                        "--",
                        x,
                        " must be lower than ",
                        length(blocks),
                        " (the maximum number of blocks), not be equal to ",
                        opt[[x]],
                        "."
                    ),
                    exit_code = 133
                )
            else if (opt[[x]] == 0)
                opt[[x]] <- length(blocks)
            else if (opt[[x]] < 0)
                stop(paste0("--",
                    x,
                    " must be positive, not be equal to ",
                    opt[[x]],
                    "."),
                    exit_code = 134)
            checkInteger(x)
        }
    }
    
    out <- lapply(seq_len(length(opt$ncomp)), function(x) {
        checkInteger("ncomp", opt$ncomp[x])
        if ((opt$ncomp[x] < 2) ||
                (opt$ncomp[x] > ncol(blocks[[x]]))) {
            stop(
                paste0(
                    "--ncomp must be comprise between 2 and ",
                    ncol(blocks[[x]]),
                    ", the number of variables of the block (currently equals to ",
                    opt$ncomp[x],
                    ")."
                ),
                exit_code = 126
            )
        }
    })
    
    checkArgSize(blocks, opt$ncomp, "ncomp")
    
    out <- sapply(c("compx", "compy"), function(x) {
        if ((opt[[x]] < 1) || (opt[[x]] > opt$ncomp[opt$block])) {
            stop(
                paste0(
                    "--",
                    x,
                    " is currently equals to ",
                    opt[[x]],
                    " and must be comprise between 1 and ",
                    opt$ncomp[opt$block],
                    " (the number of component for the selected block)."
                ),
                exit_code = 128
            )
        }
    })
    
    MSG <- "--tau must be comprise between 0 and 1 or must correspond to the character 'optimal' for automatic setting"
    if (all(opt$tau != "optimal")) {
        tryCatch({
            list_tau <- as.list(opt$tau)
            # Check value of each tau
            out <- lapply(list_tau, function(x) {
                if (((x < 0) || (x > 1)) && x != "optimal")
                    stop(paste0(MSG, " (currently equals to ", x, ")."),
                        exit_code = 129)
            })
            
            # If there is only one common tau
            if (length(list_tau) == 1)
                opt$tau <- rep(list_tau[[1]], length(blocks))
            else if (checkArgSize(blocks, list_tau, "tau"))
                opt$tau <- unlist(list_tau)
        }, warning = function(w) {
            stop(MSG, exit_code = 131)
        })
    } else {
        opt$tau <- "optimal"
    }
    
    checkC1(blocks, opt$tau, opt$type)
    
    return(opt)
}

checkInteger <- function(x, y = NULL) {
    # Test either x is an integer x : a string corresponding to a name
    # in a list opt'
    
    if (is.null(y))
        y <- opt[[x]]
    
    # Test if not a character
    tryCatch({
        as.integer(y)
        if (is.na(y)) {
            y <- opt[[x]]
            warning("")
        }
    }, warning = function(w) {
        stop(paste0("--", x, " is a character (", y, ") and must be an integer."),
            exit_code = 136)
    })
    
    # Test if not a float
    
    if (length(strsplit(as.character(y), ".", fixed = TRUE)[[1]]) > 1)
        stop(paste0("--", x, " is a float (", y, ") and must be an integer."))
}

loadLibraries <- function(librairies) {
    for (l in librairies) {
        if (!(l %in% installed.packages()[, "Package"]))
            utils::install.packages(l, repos = "http://cran.us.r-project.org")
        library(
            l,
            character.only = TRUE,
            warn.conflicts = FALSE,
            quietly = TRUE
        )
    }
}

########## Main ##########

# Get arguments : R packaging install, need an opt variable with associated
# arguments
opt <- list(
    directory = ".",
    separator = "\t",
    type = "rgcca",
    ncomp = 2,
    tau = "optimal",
    scheme = "factorial",
    init = 1,
    block = 0,
    compx = 1,
    compy = 2,
    nmark = 100,
    o1 = "individuals.pdf",
    o2 = "corcircle.pdf",
    o3 = "top_variables.pdf",
    o4 = "ave.pdf",
    o5 = "design.pdf",
    o6 = "individuals.tsv",
    o7 = "variables.tsv",
    o8 = "rgcca.result.RData",
    datasets = paste0("inst/extdata/",
        c("agriculture","industry","politic"),
        ".tsv",
        collapse = ",")
)


loadLibraries(c("RGCCA", "ggplot2", "optparse", "scales", "igraph", "ggrepel"))

tryCatch({
    opt <- checkArg(parse_args(getArgs()))
}, error = function(e) {
    if (length(grep("nextArg", e[[1]])) != 1)
        stop(e[[1]], exit_code = 140)
}, warning = function(w) {
    stop(w[[1]], exit_code = 141)
})

# Load functions
setwd(opt$directory)

for (f in list.files("R/")) {
    if ( f != "launcher.R")
        source(paste0("R/", f))
}

# Set missing parameters by default
opt$header <- !("header" %in% names(opt))
opt$superblock <- !("superblock" %in% names(opt))
# opt$bias <- !('bias' %in% names(opt))
opt$scale <- !("scale" %in% names(opt))
opt$text <- !("text" %in% names(opt))

blocks <- setBlocks(opt$datasets, opt$names, opt$separator)
blocks <- scaling(blocks, opt$scale)

# If single values for ncomp and tau, tansform it in a list
for (x in c("ncomp", "tau")) {
    if (length(grep(",", opt[[x]])) == 0 &&
            !(x == "tau" && opt[[x]] == "optimal"))
        opt[[x]] <- paste(rep(opt[[x]], length(blocks)), collapse = ",")
}

opt <- checkSuperblock(opt)
opt <- postCheckArg(opt, blocks)

if (!is.null(opt$response)) {
    opt <- setPosPar(opt, blocks, opt$response)
    blocks <- opt$blocks
}

blocks <- setSuperblock(blocks, opt$superblock, opt$type)

connection <- opt$connection
if (!is.matrix(connection))
    connection <- setConnection(blocks,
        (opt$superblock | !is.null(opt$response)),
        opt$connection,
        opt$separator)

group <- setResponse(blocks, opt$group, opt$separator, opt$header)

rgcca.out <- rgcca.analyze(
        blocks = blocks,
        connection = connection,
        tau = opt$tau,
        ncomp = opt$ncomp,
        scheme = opt$scheme,
        scale = FALSE,
        type = opt$type
    )

########## Plot ##########

# Samples common space
if (opt$ncomp[opt$block] == 1 && is.null(opt$block_y)) {
    warning("With a number of component of 1, a second block should be chosen
    to perform a samples plot")
} else {
    (
        individual_plot <- plotSamplesSpace(
                rgcca.out,
                group,
                opt$compx,
                opt$compy,
                opt$block,
                opt$text,
                opt$block_y,
                getFileName(opt$group)
            )
    )
    savePlot(opt$o1, individual_plot)
}

if (opt$ncomp[opt$block] > 1) {
    # Variables common space
    (
        corcircle <- plotVariablesSpace(
                rgcca.out,
                blocks,
                opt$compx,
                opt$compy,
                opt$superblock,
                opt$block,
                opt$text,
                n_mark = opt$nmark
            )
    )
    savePlot(opt$o2, corcircle)
}

# Fingerprint plot
top_variables <- plotFingerprint(
        rgcca.out, 
        blocks, 
        opt$compx, 
        opt$superblock, 
        opt$nmark,
        type = "cor"
    )
savePlot(opt$o3, top_variables)


if (opt$type != "pca") {
    # Average Variance Explained
    (ave <- plotAVE(rgcca.out))
    savePlot(opt$o4, ave)

    # Creates design scheme
    nodes <- getNodes(blocks, rgcca = rgcca.out)
    edges <- getEdges(connection, blocks)
    conNet <- function()plotNetwork(nodes, edges, blocks)
    savePlot(opt$o5, conNet)
}

saveInds(rgcca.out, blocks, 1, 2, opt$o6)
saveVars(rgcca.out, blocks, 1, 2, opt$o7)
save(rgcca.out, file = opt$o8)
