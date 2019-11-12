# Author: Etienne CAMENEN
# Date: 2019
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

checkNbBlocks <- function(blocks, type) {
    if (tolower(type) == "pca") {
        msg <- "Only one block is"
        exit_code <- 110
    } else{
        msg <- "Two blocks are"
        exit_code <- 111
    }

    stop(
        paste0(
            length(blocks),
            " blocks used in the analysis. ",
            msg ,
            " required for a ",
            type,
            "."
        ),
        exit_code = exit_code
    )
}

#' Define the analysis parameters
#'
#' Define the correct parameters according to the type of the analysis
#'
#' @inheritParams plotVariablesSpace
#' @param opt An OptionParser object
#' @param connection A matrix giving the connection between the blocks
#' @param tau A vector of float (or character for 'optimal' setting) giving the
#' shrinkage parameter for covariance maximization
#' @param ncomp A vector of integer giving the number of component for each
#' blocks
#' @param scheme A character giving the link function for covariance maximization
#' @param type A character giving the type of analysis
#' @param verbose A boolean displaying the warnings
#' @param quiet A boolean hidding the warnings
#' @return \item{blocks}{A list of matrix}
#' @return \item{scheme}{A character giving the link function for covariance
#' maximization}
#' @return \item{tau}{A vector of float (or character for 'optimal' setting) giving
#' the shrinkage parameter for covariance maximization}
#' @return \item{ncomp}{A vector of integer giving the number of component for each
#' blocks}
#' @return \item{connection}{matrix giving the connection between the blocks}
#' @return \item{superblock}{A boolean giving the presence (TRUE) / absence (FALSE)
#' of a superblock}
#' @export
select_type <- function(
    blocks = blocks,
    opt = NULL,
    connection = 1 - diag(length(blocks)),
    tau = rep(1, length(blocks)),
    ncomp = rep(1, length(blocks)),
    scheme = "centroid",
    superblock = TRUE,
    type  = "rgcca",
    verbose = TRUE,
    quiet = FALSE) {
    
        if (!is.null(opt)) {
            scheme <- opt$scheme
            connection <- opt$connection
            superblock <- opt$superblock
            type <- opt$type
            tau <- opt$tau
            ncomp <- opt$ncomp
            ncomp <- unlist(
                lapply(
                    strsplit(gsub(" ", "", as.character(ncomp)), ","),
                        function(x)
                            tryCatch({
                                as.double(x)
                            }, warning = function(w) {
                                stop(unique(
                                    paste0(
                                        "--ncomp is a character (",
                                        x,
                                        ") and must be an integer."
                                    )
                                ), exit_code = 136)
                            }))[[1]])
            
            l_tau <- as.list(strsplit(gsub(" ", "", as.character(tau)), ",")[[1]])

            tau <- unlist(lapply(l_tau, function(x) {
                tryCatch({
                    as.double(x)
                }, warning = function(w) {
                    "optimal"
                })
            }))
        }

        J <- length(blocks)
        MSG_SUPER <- "a superbloc is used"
        MSG_TYPE <- paste0("By using a ", toupper(type), ", ")
        warn.type.value <- warn.type.par <- warn.msg.super <- character(0)

        if (quiet)
            verbose <- FALSE

        ### SETTINGS ###

        warnParam <- function(param, x) {
            warn.type.par <<-c(warn.type.par, paste(deparse(substitute(param))))
            warn.type.value <<- c(warn.type.value, toString(x))
        }

        setTau <- function(x) {
            warnParam(tau, x)
            return(x)
        }

        setScheme <- function(x) {
            warnParam(scheme, x)
            return(x)
        }

        setConnection <- function(x) {
            warnParam(connection, paste(deparse(substitute(x))))
            return(x)
        }

        warnSuper <- function(x) {
            if (length(x) < (length(blocks))) {
                warn.msg.super <<- c(warn.msg.super, deparse(substitute(x)))
                return(c(x, x[1]))
            } else
                return(x)
        }

        setSuperbloc <- function(verbose = TRUE) {
            blocks <<- c(blocks, Superblock = list(Reduce(cbind, blocks)))
            superblock <<- TRUE
            connection <<- NULL
            ncomp <<- warnSuper(ncomp)
        }

        set2Block <- function(type) {
            if (length(blocks) != 2)
                checkNbBlocks(blocks, type)

            scheme <<- setScheme("horst")
            connection <<- setConnection(1 - diag(2))
        }

        ### CHECK TYPES ###

        if (length(grep("[sr]gcca", tolower(type))) == 1) {
            if (superblock) {
                setSuperbloc(FALSE)
                tau <- warnSuper(tau)
            } else
                superblock <- FALSE
        } else
            superblock <- FALSE

        if (length(grep("pls-?pm", tolower(type))) == 1) {
            scheme   <- setScheme("centroid")
            tau      <- setTau(rep(0, J))
            # TODO: superblock allowed in PLS-PM, whos gonna call : Arthur
        }

        else if (tolower(type) == "pca") {
            if (length(blocks) != 1)
                checkNbBlocks(blocks, type)

            scheme   <- setScheme("horst")
            tau      <- setTau(c(1, 1))
            setSuperbloc()
        }

        # 2 Blocks cases
        else if (tolower(type) %in% c("cca", "ra", "ifa", "pls")) {
            set2Block(type)

            if (tolower(type) == "cca")
                tau <- setTau(c(0, 0))

            else if (tolower(type) %in% c("ifa", "pls"))
                tau <- setTau(c(1, 1))

            else if (tolower(type) == "ra")
                tau <- setTau(c(1, 0))

        }

        # Design with 1 values everywhere
        else if (tolower(type) %in% c("sumcor",
                                    "ssqcor",
                                    "sabscor",
                                    "sumcov",
                                    "sumcov-1",
                                    "maxbet",
                                    "sabscov")) {
            connection <- setConnection(matrix(1, J, J))

            # COR models
            if (tolower(type) %in% c("sumcor", "ssqcor", "sabscor")) {
                tau <- setTau(rep(0, J))

                switch(
                    tolower(type),
                    "sumcor" = {
                        scheme <- setScheme("horst")
                    },
                    "ssqcor" = {
                        scheme <- setScheme("factorial")
                    },
                    "sabscor" = {
                        scheme <- setScheme("centroid")
                    }
                )
            }

            # COV models
            else if (tolower(type) %in% c(
                "sumcov",
                "sumcov-1",
                "maxbet",
                "ssqcov",
                "ssqcov-1",
                "maxbet-b",
                "sabscov",
                "sabscov-1"
            )) {
                tau      <- setTau(rep(1, J))

                if (tolower(type) %in% c("sumcov", "sumcov-1", "maxbet"))
                    scheme   <- setScheme("horst")

                else if (tolower(type) %in% c("ssqcov", "ssqcov-1", "maxbet-b"))
                    scheme   <- setScheme("factorial")

                else if (tolower(type) %in% c("sabscov", "sabscov-1"))
                    scheme   <- setScheme("centroid")

            }

            # Design with 1 values everywhere and 0 on the diagonal
        }

        else if (tolower(type) %in% c("sumcov-2",
                                    "maxdiff",
                                    "ssqcov",
                                    "ssqcov-1",
                                    "maxbet-b",
                                    "ssqcov-2",
                                    "maxdiff-b")) {
            connection <- setConnection(1 - diag(J))

            if (tolower(type) %in% c("sumcov-2", "maxdiff")) {
                scheme   <- setScheme("horst")
                tau      <- setTau(rep(0, J))
            }

            else if (tolower(type) %in% c("ssqcov-2", "maxdiff-b")) {
                scheme   <- setScheme("factorial")
                tau      <- setTau(rep(1, J))
            }

        }

        # Models with a superblock
        else if (tolower(type) %in% c(
            "maxvar-b",
            "gcca",
            "niles",
            "maxvar",
            "hpca",
            "maxvar-a",
            "cpca",
            "cpca-w",
            "mfa",
            "sum-pca",
            "mcoa",
            "rcon-pca",
            "ridge-gca",
            "r-maxvar"
        )) {
            setSuperbloc()

            if (tolower(type) %in% c("maxvar-b", "gcca", "niles", "maxvar")) {
                scheme   <- setScheme("factorial")
                tau      <- setTau(rep(0, J + 1))
            }

            else if (tolower(type) == "hpca") {
                scheme   <- function(x)
                    x ^ 4
                tau      <- setTau(c(rep(1, J), 0))
            }

            else if (tolower(type) %in% c(
                "maxvar-a",
                "cpca",
                "cpca-w",
                "mfa",
                "sum-pca",
                "mcoa"
            )) {
                scheme   <- setScheme("factorial")
                tau      <- setTau(c(rep(1, J), 0))
            }

            #TODO: verify these three last algo parameters

            else if (tolower(type) == "rcon-pca")
                tau <- warnSuper(tau)

            else if (tolower(type) == "ridge-gca") {
                scheme   <- setScheme("factorial")
                tau      <- setTau(c(tau[seq_len(J)], 0))
            }

            else if (tolower(type) == "r-maxvar") {
                scheme   <- setScheme("factorial")
                tau <- warnSuper(tau)
            }

        }

        else if (length(grep("[sr]gcca", tolower(type))) != 1) {
            stop(
                "Wrong type of analysis. Please select one among the following
                list: rgcca, cpca-w, gcca, hpca, maxbet-b, maxbet, maxdiff-b,
                maxdiff, maxvar-a, maxvar-b, maxvar, niles, r-maxvar, rcon-pca,
                ridge-gca, sabscor, ssqcor, ssqcor, ssqcov-1, ssqcov-2, ssqcov,
                sum-pca, sumcor, sumcov-1, sumcov-2, sumcov., sabscov, plspm.",
                exit_code = 112
            )
        }

        ### WARNINGS ###

        n = length(warn.type.par)
        if (verbose & n > 0) {
            setPlural = function(x = warn.type.par,
                                y = warn.type.value,
                                sep = " and ") {
                warn.type.par <<- paste0(x, collapse = sep)
                warn.type.value <<- paste0(y, collapse = sep)
            }

            if (n > 1) {
                grammar = "s were respectively"
                if (n == 2)
                    setPlural()
                else{
                    warn.type = c(warn.type.par[n], warn.type.value[n])
                    setPlural(warn.type.par[-n], warn.type.value[-n], ", ")
                    setPlural(c(warn.type.par, warn.type[1]),
                            c(warn.type.value, warn.type[2]))
                }
            } else
                grammar <- " was"

            msg <- paste0(warn.type.par,
                        " parameter",
                        grammar,
                        " set to ",
                        warn.type.value)

            if (superblock & tolower(type) != "pca")
                msg <- paste0(msg, " and ", MSG_SUPER)

            warning(paste0(MSG_TYPE, msg , "."))
        }

        if (verbose & superblock) {
            if (n < 0)
                paste0(MSG_SUPER, MSG_SUPER)
        }

        if (!quiet & length(warn.msg.super) > 0) {
            if (length(warn.msg.super) > 1) {
                warn.msg.super <- paste(warn.msg.super, collapse = " and ")
                grammar <- "were those"
            } else
                grammar <- "was the one"

            # warning(paste0("By using a superblock, ", warn.msg.super,
            #    " of the superblock ", grammar," of the first block."))
        }

        opt$blocks <- blocks
        opt$scheme <- scheme
        opt$tau <- tau
        opt$ncomp <- ncomp
        opt$connection <- connection
        opt$superblock <- superblock

        return(opt)
    }

#' Performs a r/sgcca
#'
#' Performs a r/sgcca with predefined parameters
#' @inheritParams select_type
#' @param scale A boolean scaling the blocks
#' @param init A character among "svd" (Singular Value Decompostion) or "random"
#' for alorithm initialization
#' @param bias A boolean for a biased variance estimator
#' @param type A character giving the type of analysis
#' @param verbose A boolean to display the progress of the analysis
#' @return A RGCCA object
#' @examples 
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.analyze(blocks)
#' @export
rgcca.analyze <- function(
    blocks,
    connection = 1 - diag(length(blocks)),
    tau = rep(1, length(blocks)),
    ncomp = rep(2, length(blocks)),
    scheme = "factorial",
    scale = TRUE,
    init = "svd",
    bias = TRUE,
    type = "rgcca",
    verbose = TRUE) {

    WARN <- FALSE

    for (i in seq_len(length(blocks))) {
        if (ncol(blocks[[i]]) > 1000) {
            # if( (type <-<- "sgcca" && tau > 0.3) || type !<- "sgcca" )
            WARN <- TRUE
        }
    }

    if (WARN & verbose)
        message("RGCCA in progress ...")

    if (tolower(type) == "sgcca") {
        func <- sgcca
        par <- "c1"
    } else{
        func <- rgcca
        par <- "tau"
    }

    func.complete <- quote(
        func(
            A = blocks,
            C = connection,
            scheme = scheme,
            ncomp = ncomp,
            scale = scale,
            verbose = FALSE,
            init = init,
            bias = bias
        )
    )
    func.complete[[par]] <- tau

    func.res <- eval(as.call(func.complete))
    names(func.res$a) <- names(blocks)

    return(func.res)
}

#' Compute bootstrap (internal)
#' 
#' Internal function for computing boostrap of RGCCA
#' 
#' @inheritParams rgcca.analyze
#' @inheritParams plotVariablesSpace
#' @return A list of RGCCA bootstrap weights
#' @examples 
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' bootstrap_k(blocks, rgcca.res, FALSE)
#' @export
bootstrap_k <- function(
    blocks,
    rgcca,
    scale = TRUE,
    init = "svd",
    bias = TRUE) {

    # Shuffle rows
    id_boot <- sample(NROW(blocks[[1]]), replace = TRUE)

    if (isTRUE(scale))
            boot_blocks <- lapply(blocks, function(x)
            scale(
                x[id_boot, ],
                center = attr(blocks, "scaled:center"),
                scale = attr(blocks, "scaled:scale")
            ) / sqrt(ncol(x)))
    else
        boot_blocks <- lapply(blocks, function(x)
            scale2(x[id_boot, ], scale = FALSE))

    boot_blocks <- removeColumnSdNull(boot_blocks)
    
    if (is(rgcca, "sgcca"))
        tau <- rgcca$c1
    else
        tau <- rgcca$tau

    # Get boostraped weights
    w <- rgcca.analyze(
        boot_blocks,
        rgcca$C,
        tau = tau,
        ncomp = rgcca$ncomp,
        scheme = rgcca$scheme,
        scale = FALSE,
        init = init,
        bias = bias,
        type = class(rgcca),
        verbose = FALSE
    )$a

    # Add removed variables
    missing_var <- lapply(seq_len(length(blocks)), function(x)
        setdiff(colnames(blocks[[x]]), rownames(w[[x]])))
    missing_tab <- lapply(seq_len(length(missing_var)),
                        function(x)
                            matrix(
                                0,
                                length(missing_var[[x]]),
                                rgcca$ncomp[x],
                                dimnames = list(missing_var[[x]], seq_len(rgcca$ncomp[x]))
                            ))

    # bug mapply with pca
    w <- lapply(seq_len(length(w)), function(x)
        rbind(w[[x]], missing_tab[[x]]))
    w <- lapply(seq_len(length(w)), function(x)
        w[[x]][colnames(blocks[[x]]), ])

    names(w) <- names(blocks)
    return(w)
}

#' Compute bootstrap
#' 
#' Computing boostrap of RGCCA
#' 
#' @inheritParams rgcca.analyze
#' @inheritParams plotVariablesSpace
#' @param n_boot A integer for the number of boostrap
#' @param nb_cores An integer for the number of cores used in parallelization
#' @return A list of RGCCA bootstrap weights
#' @examples 
#' library(RGCCA) 
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' bootstrap(blocks, rgcca.res, 2, FALSE)
#' @export
bootstrap <- function(
    blocks,
    rgcca,
    n_boot = 5,
    scale = TRUE,
    init = "svd",
    bias = TRUE,
    nb_cores = parallel::detectCores() - 1) {

    if (nb_cores == 0)
        nb_cores <- 1

    if (any(unlist(lapply(blocks, ncol) > 1000)))
        verbose <- TRUE

    w1 <- rgcca$a

    cat("Bootstrap in progress...")

    W <- parallel::mclapply(seq_len(n_boot), function(x) {

        w <- bootstrap_k(blocks,
                        rgcca,
                        scale,
                        init,
                        bias)

        # Test on the sign of the correlation
        for (k in seq_len(length(blocks))) {
            for (j in seq_len(ncol(w[[k]]))) {
                if (cor(w1[[k]][, j], w[[k]][, j]) < 0)
                    w[[k]][, j] <- -1 * w[[k]][, j]
            }
        }

        return(w)

    }, mc.cores = nb_cores)

    cat("OK", append = TRUE)

    return(W)
}

#' Extract a bootstrap
#' 
#' Extract statistical information from a bootstrap
#'
#' @inheritParams bootstrap
#' @inheritParams plotHistogram
#' @inheritParams plotVariablesSpace
#' @param W A list of list weights (one per bootstrap per blocks)
#' @param comp An integer giving the index of the analysis components
#' @return A matrix containing the means, 95% intervals, bootstrap ratio and p-values
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' boot = bootstrap(blocks, rgcca.res, 2, FALSE)
#' getBootstrap(rgcca.res, boot)
#' @export
getBootstrap <- function(
    rgcca,
    W,
    comp = 1,
    i_block = NULL,
    collapse = TRUE,
    nb_cores = parallel::detectCores() - 1) {

    if (nb_cores == 0)
        nb_cores <- 1

    if (is.null(i_block))
        i_block <- length(W[[1]])

    if (comp > min(rgcca$ncomp))
        stop("Selected dimension was not associated to every blocks",
            exit_code = 113)

    cat("Binding in progress...")
    
    mean <- weight <- sd <- occ <- list()

    if (collapse)
        J <- seq(length(rgcca$a))
    else
        J <- i_block

    for (i in J) {

        W_bind <- parallel::mclapply(
            W,
            function(x) x[[i]][, comp],
            mc.cores = nb_cores
        )
        
        weight[[i]] <- rgcca$a[[i]][, comp]
        W_select <- matrix(unlist(W_bind), nrow = length(W_bind), ncol = length(W_bind[[1]]), byrow = TRUE)
        colnames(W_select) <- names(weight[[i]])
        rm(W_bind); gc()

        n <- seq(ncol(W_select))

        if (is(rgcca, "sgcca")) {

            occ[[i]] <- unlist(parallel::mclapply(n,
                function(x) sum(W_select[,x] != 0) / length(W_select[, x]),
                mc.cores = nb_cores
            ))

        }
            
        mean[[i]] <- unlist(parallel::mclapply(n,
            function(x) mean(W_select[,x]),
            mc.cores = nb_cores
        ))
        sd[[i]] <- unlist(parallel::mclapply(n,
            function(x) sd(W_select[,x]),
            mc.cores = nb_cores
        ))

        rm(W_select); gc()
    }

    rm(W); gc()
    
    occ <- unlist(occ)
    mean <- unlist(mean)
    weight <- unlist(weight)
    sd <- unlist(sd)

    cat("OK", append = TRUE)

    p.vals <- pnorm(0, mean = abs(mean), sd = sd)
    tail <- qnorm(1 - .05 / 2)

    df <- data.frame(
        mean = mean,
        rgcca = weight,
        intneg = mean - tail * sd,
        intpos = mean + tail * sd,
        br = abs(mean) / sd,
        p.vals,
        BH = p.adjust(p.vals, method = "BH")
    )

    if (is(rgcca, "sgcca")) {
        index <- 8
        df$occ <- occ
    }else{
        index <- 5
        df$sign <- rep("", nrow(df))

        for (i in seq(nrow(df)))
            if (df$intneg[i]/df$intpos[i] > 0)
                df$sign[i] <- "*"
    }

    if (collapse)
        df$color <- as.factor(getBlocsVariables(rgcca$a, collapse = collapse))
    
    zero_var <- which(df[, 1] == 0)
    if (length(zero_var) != 0)
        df <- df[-zero_var, ]
    
    data.frame(getRankedValues(df, index, allCol = TRUE), order = nrow(df):1)
}

#' Plot a bootstrap
#' 
#' Plot the top variables from a bootstrap
#'
#' @inheritParams plotHistogram
#' @inheritParams plotVariablesSpace
#' @param show.boot A boolean to show the bootstrap mean and sd on the graphic
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' boot = bootstrap(blocks, rgcca.res, 2, FALSE)
#' selected.var = getBootstrap(rgcca.res, boot)
#' plotBootstrap(selected.var, rgcca.res)
#' @export
plotBootstrap <- function(
    df,
    rgcca,
    superblock = TRUE,
    show.boot = TRUE,
    n_mark = 30) {

    color <- intneg <- intpos <- NULL
    J <- names(rgcca$a)

    if (nrow(df) > n_mark)
        df <- df[seq_len(n_mark), ]

    if (superblock) {
        color2 <- factor(df$color)
        levels(color2) <- colorGroup(color2)
        p <- ggplot(df, aes(order, mean, fill = color))
    } else{
        color2 <- "black"
        p <- ggplot(df, aes(order, mean, fill = abs(mean)))
    }

    p <- plotHistogram(p, df, "Variable mean", as.character(color2))
    
    if (show.boot) {
    p <- p +
        geom_line(aes(x = order, y = mean), inherit.aes = FALSE, lwd = 0.7) +
        geom_point(aes(x = order, y = mean), inherit.aes = FALSE, size = 1.5)

        if (is(rgcca, "rgcca" ))
            p <- p +
                geom_errorbar(aes(ymin = intneg, ymax = intpos))
    }

    if (superblock)
        col <- J
    else
        col <- J[-length(J)]
        
    if (superblock) {
        matched <- match(rev(unique(df$color)), col)
        p <- orderColorPerBlocs(rgcca$a, p, matched, superblock)
    }

    return(p)
}

#' Plot a bootstrap in 2D
#' 
#' Biplot of the top variables from a SGCCA bootstrap with the number of 
#' non-zero occurences in x-axis and the boot-ratio (mean/sd) in y-axis. 
#' Negative weights are colored in red and the positive ones are in green.
#'
#' @param b A matrix of boostrap
#' @param x A character for the column to plot in x-axis
#' @param y A character for the column to plot in y-axis
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks, type = "sgcca")
#' boot = bootstrap(blocks, rgcca.res, 2, FALSE)
#' selected.var = getBootstrap(rgcca.res, boot)
#' plotBootstrap2D(selected.var)
#' @export
plotBootstrap2D <- function(b, x = "br", y = "occ"){

    axis <- function(margin){
        element_text(
            face = "italic",
            size = AXIS_TITLE_CEX * 0.75,
            margin = margin
        )
    }

    ggplot(
        b, 
        aes(
            x = abs(b[, x]),
            y = b[, y],
            label = row.names(b),
            color = as.factor(mean > 0)
        )
    ) + 
    geom_text(
        size = PCH_TEXT_CEX * 0.75
    ) + 
    labs(
        y = "Non-zero occurences",
        x = "Bootstrap-ratio", 
        title = "Occurences selection\nby bootstrap"
    ) + 
    theme_classic()  +
    theme_perso() +
    theme(
        legend.position = "none",
        axis.title.y = axis(margin(0, 20, 0, 0)),
        axis.title.x = axis(margin(20, 0, 0, 0))
    ) +
    scale_color_manual(values = colorGroup(seq(2)))
}

#' Plot a bootstrap in 1D
#' 
#' Histogram of the best variables of an RGCCA bootstrap with, on the x-axis, 
#' the number of non-zero occurrences (SGCCA) or the bootstrap-ratio 
#' (mean/sd; RCCA). The bars are colored according to the average weight of 
#' the boostrap  (according to an ascending gradient from red to blue)
#'
#' @param b A matrix of boostrap
#' @param x A character for the column used in the plot
#' @param y A character for the column to color the bars
#' @param n An integer giving the number maximum of top variables
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' boot = bootstrap(blocks, rgcca.res, 2, FALSE)
#' selected.var = getBootstrap(rgcca.res, boot)
#' plotBootstrap1D(selected.var)
#' @export
plotBootstrap1D <- function(b, x = "occ", y = "mean", n = 50){

    if (!("occ" %in% colnames(b))) {
        title <- "Bootstrap ratio"
        x <- "br"
    }else
        title <- "Occurences selection\nby bootstrap"

    b <- head(b, n)
    p <- ggplot(b, 
        aes(
            x = order,
            y = b[, x],
            fill = b[, y]
        )
    )

    plotHistogram(
        p,
        b,
        title,
        "black",
        low_col = colorGroup(seq(3))[1],
        mid_col = "white", 
        high_col = colorGroup(seq(3))[3]
    ) +
    labs(fill = "Mean weights")
}

scaling = function(blocks,
                    scale = TRUE,
                    bias = TRUE) {
    if (scale) {
        lapply(blocks, function(x)
            scale2(x, bias = bias) / sqrt(ncol(x)))
    }else{
        lapply(blocks, function(x)
            scale2(x, scale = FALSE))
    }
}

checkC1 <- function(blocks, tau, type) {
    # c1 : A vector of integer giving the spasity parameter for SGCCA (c1)
    # Stop the program if at least one c1 parameter is not in the required interval

    if (tolower(type) == "sgcca") {
        #the minimum value avalaible
        min_c1 <- lapply(blocks, function(x) 1 / sqrt(ncol(x)))

        # Check c1 varying between 1/sqrt(pj) and 1
        mapply(function(x, y) {
            if (x < y | x > 1)
                stop(
                    paste0(
                        "Sparsity parameter is equals to ",
                        x,
                        ". For SGCCA, it must be comprise between 1/sqrt(number_column) (i.e., ",
                        toString(unlist(
                            lapply(min_c1, function(x)
                                ceiling(x * 100) / 100)
                        ))
                        ,
                        ") and 1."
                    ),
                    exit_code = 132
                )
        }, tau, min_c1)
    }
}
