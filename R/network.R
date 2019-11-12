# Author: Etienne CAMENEN
# Date: 2019
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

#' Creates the nodes for a design matrix
#' 
#' @inheritParams plotVariablesSpace
#' @inheritParams select_type
#' @return A dataframe with blocks in rows and the number of variables, of rows 
#' and tau or c1 in columns
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' getNodes(blocks, rgcca = rgcca.res)
#' @export
getNodes <- function(blocks, tau = NULL, rgcca = NULL) {

    if (!is.null(rgcca) & is(rgcca, "sgcca")) {
        par.rgcca <- "c1"
        par.name <- "sparsity"
    } else
        par.rgcca <- par.name <- "tau"

    if (any(tau == "optimal")) {
        if (!is.null(rgcca))
            tau <- unlist(lapply(seq_len(ncol(rgcca[[par.rgcca]])), function(x)
                    Reduce(paste, round(rgcca[[par.rgcca]][, x], 2))))
        else
            tau <- rep(NA, length(blocks))
    }

    if (is.null(tau)) {
        if (is.matrix(rgcca[[par.rgcca]]))
            tau <-  unlist(lapply(seq_len(ncol(rgcca[[par.rgcca]])), function(x)
                    Reduce(paste, round(rgcca[[par.rgcca]][, x], 2))))
        else
            tau <- rgcca[[par.rgcca]]
    }

    nrow <- unlist(lapply(blocks, function(x)
            ifelse(
                is.null(attributes(x)$nrow),
                nrow(blocks[[1]]),
                attributes(x)$nrow
            )))

    values <- list(names(blocks), unlist(lapply(blocks, NCOL)), nrow, tau)
    nodes <- as.data.frame(matrix(unlist(values), length(blocks), length(values)))
    colnames(nodes) <- c("id", "P", "nrow", par.name)

    return(nodes)
}

#' Creates the edges for a design matrix
#' 
#' @inheritParams select_type
#' @return A dataframe with tuples of connected blocks
#' @examples
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' getEdges(rgcca.res$C, blocks)
#' @export
getEdges <- function(connection, blocks) {
    J <- NCOL(connection)

    edges <- list()

    k <- 0
    for (j in seq_len(J)) {
        for (i in seq_len(J)) {
            if (i > k && connection[i, j] > 0)
                edges[[length(edges) + 1]] <-
                    c(names(blocks)[j], names(blocks)[i], connection[i, j])
        }
        k <- k + 1
    }

    edges <- as.data.frame(t(matrix(unlist(edges), 3, length(edges))))
    colnames(edges) <- c("from", "to", "weight")
    edges[, 3] <- as.numeric(edges[, 3])

    return(edges)
}

colorNodes <- function(nodes) {
    unlist(lapply(as.list(1 - nodes$P / max(nodes$P)), function(x)
        rgb(colorRamp( c("coral3", "khaki2"))(x) / 255)))
}

#' Plot the connection between blocks
#' 
#' @inheritParams select_type
#' @param nodes A dataframe containing metadata for each blocks
#' @param edges A dataframe of connection between blocks
#' @return A dataframe with tuples of connected blocks
#' @examples
#' library(igraph)
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' e <- getEdges(rgcca.res$C, blocks)
#' n <- getNodes(blocks, rgcca = rgcca.res)
#' plotNetwork(n, e, blocks)
#' @export
plotNetwork <- function(nodes, edges, blocks) {
    # Avoid random
    set.seed(1)
    V <- E <- NULL

    par <- ifelse("sparsity" %in% names(nodes), "sparsity", "tau")

    net <- graph_from_data_frame(d = edges,
            vertices = nodes,
            directed = FALSE)

    if (all(is.na(nodes[, par]))) {
        nodes[, par] <- rep("optimal", length(blocks))
        V(net)$tau <- rep(1, length(blocks))
    }

    V(net)$color <- "khaki2"
    V(net)$label <-
        paste(nodes$id,
            "\nP =",
            nodes$P,
            paste0("\n", par, " ="),
            nodes[,par],
            "\nN =",
            nodes$nrow,
            sep = " ")
    V(net)$shape <- "square"
    E(net)$width <- E(net)$weight * 2

    plot(
        net,
        cex.main = 5,
        edge.color = "gray70",
        edge.lty = 2,
        vertex.frame.color = "gray50",
        vertex.label.color = "black",
        vertex.label.dist = 6,
        vertex.label.degree = 1.5,
        vertex.size = 23,
        main = paste0("Common rows between blocks : ", nrow(blocks[[1]]))
    )

}

#' Plot the connection between blocks (dynamic plot)
#' 
#' @inheritParams select_type
#' @inheritParams plotNetwork
#' @return A dataframe with tuples of connected blocks
#' @examples
#' library(visNetwork)
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' e <- getEdges(rgcca.res$C, blocks)
#' n <- getNodes(blocks, rgcca = rgcca.res)
#' plotNetwork2(n, e, blocks)
#' @export
plotNetwork2 <- function(nodes, edges, blocks) {
    par <- ifelse("sparsity" %in% names(nodes), "sparsity", "tau")

    if (all(is.na(nodes[, par])))
        nodes[, par] <- rep("optimal", length(blocks))

    nodes$title <- nodes$id
    nodes$label <- paste(nodes$id,
            "\nP =",
            nodes$P,
            paste0("\n", par, " ="),
            nodes[,par],
            "\nN =",
            nodes$nrow,
            sep = " ")

    edges$width <- edges$weight * 2
    nodes$color.background <- rep("#eee685", length(blocks))

    visNetwork(
        nodes,
        edges,
        main = list(
            text = paste0("Common rows between blocks : ",
                        nrow(blocks[[1]])),
            style = "font-family:sans;font-weight:bold;font-size:28px;text-align:center;"
        )
    ) %>%
        visNodes(
            borderWidth = 2,
            shape = "square",
            shadow = TRUE,
            color = list(
                border = "gray",
                highlight = list(background = "black", border = "darkred")
            )
        ) %>% visEdges(
            smooth = FALSE,
            shadow = TRUE,
            dashes = TRUE,
            color = list(color = "gray", highlight = "darkred")
        )

}
