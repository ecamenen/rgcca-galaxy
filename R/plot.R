# Author: Etienne CAMENEN
# Date: 2019
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

CEX <- 1
AXIS_TITLE_CEX <- 19 * CEX
SUBTITLE_CEX <- 16 * CEX
AXIS_TEXT_CEX <- 10 * CEX
PCH_TEXT_CEX <- 3 * CEX

# X- Y axis format for plotly objets : no axis, no ticks
ax <- list(linecolor = "white",
        ticks = "",
        titlefont = list(size = 23))
ax2 <- list(linecolor = "white",
        tickfont = list(size = 10, color = "grey"))

# Dynamic visualization of the outputs
# f: ggplot2 function
# ax: list object containing attributes of xaxis / yaxis parameter in plotly
# (https://plot.ly/javascript/reference/, xaxis/yaxis)
# text: axis information to print (among y, x, x+y, text)
#  (https://plot.ly/javascript/reference/, hoverinfo)
# dynamicTicks: a bolean giving the generation of axis tick labels
# (otherwhise samplesPlot which do not have traces could not be convereted
# in ggplotly)
# return a plotly object
dynamicPlot <- function(f,
    ax,
    text = "name+x+y",
    legend = TRUE,
    dynamicTicks = FALSE) {
    
        # Convert a ggplot into a plotly object add a layout with predefined
        # formats for
        # x- and y- axis set the style to show onMouseOver text
        p <- plotly_build(
                ggplotly(f, dynamicTicks = dynamicTicks) %>%
                layout(
                    xaxis = ax,
                    yaxis = ax,
                    annotations = list(showarrow = FALSE, text = "")
                ) %>% style(hoverinfo = text)
            )
        
        if (legend) {
            
            # set on the top the position of the legend title
            p$x$layout$annotations[[1]]$yanchor <- "top"
            # Deals with a too short name of modalities
            p$x$layout$margin$r <- nchar(p$x$layout$annotations[[1]]$text) * 13
            p$x$layout$margin$t <- 100
            # for shiny corcircle, if text = TRUE, two legends will appear.
            # only the first one will be selected
            title <- unlist(strsplit(p$x$layout$annotations[[1]]$text, "<br />"))[1]

            # to prevent print a 'NA' when there is no legend in plot
            if (is.na(title))
                title <- ""

            # set the font for this title
            p$x$layout$annotations[[1]]$text <- paste0("<i><b>", title, "</b></i>")
            # Sys.info()[['sysname']]

            if (!is.null(f$labels$subtitle)) {
                if (packageVersion("plotly") < 4.9)
                    p$x$layout$title <- paste0(
                            p$x$layout$title,
                            "<br><b>",
                            "c",
                            substring(f$labels$subtitle, 2),
                            "</b>"
                        )
                else
                    p$x$layout$title$text <- paste0(
                        p$x$layout$title$text,
                        "<br><b>",
                        "c",
                        substring(f$labels$subtitle, 2),
                        "</b>"
                    )
            }
        }
        
        if (ncol(f$data) == 3)
            p$sample_names <- lapply(
                levels(as.factor(f$data[, 3])), 
                function(x) row.names(subset(f$data, f$data[, 3] == x)))
        else
            p$sample_names <- list(row.names(f$data))
        return(p)
    }


dynamicPlotBoot <- function(p) {
    
    p <- dynamicPlot(p, ax2, "text")
    n <- length(p$x$data)
    m <- unlist(lapply(p$x$data, function(x) !is.null(x$orientation)))
    j <- length(m[m])
    
    for (i in seq_len(j)) {
        p$x$data[[i]]$text <- paste(
            round(p$x$data[[i]]$x, 3),
            "+/-",
            round(p$x$data[[n]]$error_x$array[j],3))

        j <- j - 1
    }

    # Remove the onMouseOver for the last trace
    changeText(p) %>% style(
        error_x = list(
            array = p$x$data[[n]]$error_x$array,
            color = "gray"
        ),
        hoverinfo = "none",
        traces = n
    )
}

# p: a ggplot function
# hovertext : a boolean for the use of hovertext (if TRUE) as the attribute 
# to parse the onMouseOver text ("text" attribute, if FALSE)
changeHovertext <- function(p, hovertext = TRUE) {
    
    attr <- ifelse(hovertext, "hovertext", "text")
    # identify the order / id of the traces which corresponds to x- and y-axis 
    # (should be before the splitting function)
    traces <- which(lapply(
        p$x$data,
        function(x) length(grep("intercept", x$text)) == 1) == TRUE)
    
    # length of groups of points without traces and circle points
    n <- which(sapply(p$x$data, function(x)
        match("xintercept: 0", x$text) == 1)) - 1

    for (i in seq_len(n)) {
        # For each lines of each group in the legend
        for (j in seq_len(length(p$x$data[[i]][attr][[1]]))) {
            
            # Distinguish each duplicate by splitting with "<br>" and separe
            # them in key/value by splitting with ": " (like a dictionnary)
            l_text <- lapply(
                as.list(strsplit(p$x$data[[i]][attr][[1]][j], "<br />")[[1]]),
                function(x) strsplit(x, ": ")[[1]])
            
            # keep only the (x, y) coordinates with the key df[, ]
            # and the response if exists
            l_text = unlist(lapply(l_text, function(x, y) {
                if (x[1] %in% paste0("df[, ", c(1, 2), "]"))
                    round(as.numeric(x[2]), 3)
                else if (x[1] == "resp")
                    x[2]
            }))

            # do not print names because text = FALSE plots have'nt names
            name = ifelse(hovertext,
                        paste0("name: ", p$x$data[[i]]$text[j], "<br />"),
                        "")
            # Overwrite the onMouseOver text with the (x, y) coordinates and
            # the response if exists
            p$x$data[[i]][attr][[1]][j] <-
                paste0(name,"x: ",l_text[1], "<br />y: ", l_text[2],
                            ifelse(
                                length(l_text) == 3,
                                paste0("<br />response: ", l_text[3]) ,
                                ""
                            ))

        }
    }
    # Remove the x- and y- axis onOverMouse
    (style(p, hoverinfo = "none", traces = traces))
}

changeText <- function(p) {
    for (i in seq_len(length(p$x$data)))
        p$x$data[[i]]$text <- sub(
            "order: .*<br />df\\[, 1\\]: (.*)<.*",
            "\\1\\",
            p$x$data[[i]]$text)
    return(p)
}

#' Print the variance of a component
#'
#' Prints the percent of explained variance for a component of a block 
#' (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @inheritParams plotSamplesSpace
#' @param n An integer giving the index of the analysis component
#' @param i An integer giving the index of a list of blocks
#' @param outer A boolean for ave plot case
#' @return A string for the variance on the component
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' AVE = list(c(0.6, 0.5), c(0.7, 0.45))
#' rgcca.res = list(AVE = list(AVE_X = AVE))
#' # For the superblock (or the last block)
#' printAxis(rgcca.res, 1)
#' # "Axis 1 (70%)"
#' # For the first block
#' printAxis(rgcca.res, 2, 1)
#' # "Axis 2 (50%)"
#' @export
printAxis <- function(rgcca, n = NULL, i = NULL, outer = FALSE) {
    
    # by default, take the last block
    if (is.null(i))
        i <- length(rgcca$AVE$AVE_X)

    nvar <- sum(rgcca$a[[i]][, n] != 0)
    if (!is(rgcca, "sgcca") | nvar == length(rgcca$a[[i]][, n]))
        varText <- ""
    else
        varText <- paste0(nvar, " variables, ")
    
    ave <- quote(paste0(round(AVE[n] * 100, 1), "%"))
    if (isTRUE(outer)) {
        AVE <- rgcca$AVE$AVE_outer
        n <- c(1, 2)
        paste0("First outer comp. : ", paste(eval(ave), collapse = " & "))
    } else {
        AVE <- rgcca$AVE$AVE_X[[i]]
        paste0("Comp. ", n, " (", varText, eval(ave), ")")
    }
}

# Default theme for ggplot
theme_perso <- function() {
    theme(
        legend.text = element_text(size = 13 * CEX),
        legend.title = element_text(face = "bold.italic", size = SUBTITLE_CEX ),
        plot.title = element_text(
            size = 25 * CEX,
            face = "bold",
            hjust = 0.5,
            margin = margin(0, 0, 20, 0)
        )
    )
}

#' Groups of color
#' 
#' Returns a color vector of equal size to the input vector
#' @param x A vector
#' @return A color vector of equal size to the input vector
#' @examples
#' colorGroup(seq(10))
#' @export
colorGroup <- function(x) {
    palette <-
        rep(
            c(
                "#cd5b45",
                "#71ad65",
                "#3c78b4",
                "#ffc600",
                "#b448af",
                "#9d9d9d",
                "#abcaef",
                "#4a6f43",
                "#f0e500",
                "#efb8f0",
                "black",
                "#d6d6d6"
            ),
            10
        )
    palette[0:length(levels(as.factor(x)))]
}

#' Plot the two components of a RGCCA
#'
#' Plot the two components of a RGCCA
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param resp A vector of characters corresponding either to a qualitative
#' variable with levels or a continuous variable
#' @param comp_x An integer giving the index of the analysis component used
#' for the x-axis
#' @param comp_y An integer giving the index of the analysis component used
#' for the y-axis
#' @param i_block An integer giving the index of a list of blocks
#' @param text A bolean to represent the points with their row names (TRUE)
#' or with circles (FALSE)
#' @param i_block_y An integer giving the index of a list of blocks (another
#' one, different from the one used in i_block)
#' @param reponse_name A character giving the legend title
#' @param no_Overlap A boolean to avoid overlap in plotted text
#' @param predicted A list containing as  2nd element a matrix of predicted components 
#' @examples
#' coord = lapply(seq_len(3),
#'    function(x) matrix(runif(15 * 2, min = -1), 15, 2))
#' AVE_X = lapply(seq_len(3), function(x) runif(2))
#' for (i in 1:length(coord))
#' row.names(coord[[i]]) = seq(15)
#' rgcca.res = list(Y = coord, AVE = list(AVE_X = AVE_X))
#' # Using a superblock
#' resp = as.matrix(rep(LETTERS[seq_len(3)], each = 5))
#' row.names(resp) = seq(15)
#' plotSamplesSpace(rgcca.res, resp)
#' # Using the first block
#' resp = as.matrix(runif(15, min=-15, max = 15))
#' row.names(resp) = seq(15)
#' plotSamplesSpace(rgcca.res, resp, 1, 2, 1)
#' @export
plotSamplesSpace <- function(
    rgcca,
    resp,
    comp_x = 1,
    comp_y = 2,
    i_block = length(rgcca$Y),
    text = TRUE,
    i_block_y = i_block,
    reponse_name = "Response",
    no_Overlap = FALSE,
    predicted = NULL) {

    if (is.null(i_block_y))
        i_block_y <- i_block

    df <- getComponents(
        rgcca = rgcca,
        resp = resp,
        comp_x = comp_x,
        comp_y = comp_y,
        i_block = i_block,
        i_block_y = i_block_y,
        predicted = predicted
    )

    if (nrow(df) > 100)
        PCH_TEXT_CEX <- 2

    if (!is.null(predicted))
            p <- ggplot(df, aes(df[, 1], df[, 2], color = df$resp))

    else if (length(unique(as.matrix(df$resp))) > 5 && 
            !unique(isCharacter(as.vector(df$resp))) ) {

        p <- ggplot(df, aes(df[, 1], df[, 2], color = df$resp))

    }else
        p <- NULL


    p <- plotSpace(rgcca,
                    df,
                    "Sample",
                    df$resp,
                    reponse_name,
                    comp_x,
                    comp_y,
                    i_block,
                    p,
                    text,
                    i_block_y,
                    no_Overlap = no_Overlap)

    # remove legend if missing
    if (length(unique(df$resp)) == 1)
        p + theme(legend.position = "none")
    else
        p
}

#' Get the blocs of each variables
#'
#' Get a vector of block names for each corresponding variable. The last block 
#' is considered as the superblock and ignored.
#'
#' @param df A list of matrix where their names are those of the blocks and the 
#' superblock and their rows are named after their variables
#' @param collapse A boolean to combine the variables of each blocks as result
#' @return A vector of character giving block names for each corresponding 
#' variable.
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' rgcca.res = list(a = rep(NA, 4))
#' names(rgcca.res$a) = LETTERS[seq_len(4)]
#' getBlocsVariables(rgcca.res)
#' # a, b, c
#' @export
getBlocsVariables <- function(df, collapse = FALSE) {
    
    if (!collapse)
        bl.names <- names(df)[-length(df)]
    else
        bl.names <- names(df)

    res <- rep(
        bl.names,
        sapply(
            df[seq_len(length(df) - as.integer(!collapse))],
            function(x) nrow(as.matrix(x))
        )
    )
    
    names(res) <- unlist(lapply(df[bl.names], row.names))

    return(res)
}

#' Plot of variables space
#'
#' Correlation circle highlighting the contribution of each variables to the
#' construction of the RGCCA components
#' @inheritParams plotSamplesSpace
#' @param blocks A list of matrix
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE)
#' of a superblock
#' @param removeVariable A bolean to keep only the 100 variables of each
#' component with the biggest correlation#'
#' @param n_mark An integer giving the number of top variables to select
#' @param collapse A boolean to combine the variables of each blocks as result
#' @examples
#' setMatrix = function(nrow, ncol, iter = 3) lapply(seq_len(iter),
#'     function(x) matrix(runif(nrow * ncol), nrow, ncol))
#' blocks = setMatrix(10, 5)
#' blocks[[4]] = Reduce(cbind, blocks)
#' for (i in seq_len(4))
#'     colnames(blocks[[i]]) = paste0( LETTERS[i],
#'     as.character(seq_len(NCOL(blocks[[i]]))))
#' coord = setMatrix(10, 2, 4)
#' a = setMatrix(5, 2)
#' a[[4]] = matrix(runif(15 * 2), 15, 2)
#' AVE_X = lapply(seq_len(4), function(x) runif(2))
#' rgcca.res = list(Y = coord, a = a, AVE = list(AVE_X = AVE_X))
#' names(rgcca.res$a) = LETTERS[seq_len(4)]
#' # Using a superblock
#' plotVariablesSpace(rgcca.res, blocks, 1, 2, TRUE)
#' # Using the first block
#' plotVariablesSpace(rgcca.res, blocks, 1, 2, FALSE, 1)
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' # Without superblock but with the of all variables to the first block
#' plotVariablesSpace(rgcca.res, blocks, collapse = TRUE)
#' @export
plotVariablesSpace <- function(
    rgcca,
    blocks,
    comp_x = 1,
    comp_y = 2,
    superblock = TRUE,
    i_block = length(blocks),
    text = TRUE,
    removeVariable = TRUE,
    n_mark = 100,
    collapse = FALSE,
    no_Overlap = FALSE) {

    y <- NULL

    df <- getVariablesIndexes(
        rgcca = rgcca,
        blocks = blocks,
        comp_x = comp_x,
        comp_y = comp_y,
        i_block = i_block,
        type = "cor",
        superblock = superblock,
        n_mark = n_mark,
        collapse = collapse,
        removeVariable = removeVariable
    )

    circleFun <- function(center = c(0, 0), diameter = 2, npoints = 100) {
        r <- diameter / 2
        tt <- seq(0, 2 * pi, length.out = npoints)
        xx <- center[1] + r * cos(tt)
        yy <- center[2] + r * sin(tt)
        return(data.frame(x = xx, y = yy))
    }
    
    if (superblock)
        collapse <- TRUE

    p <- plotSpace(
            rgcca,
            df,
            "Variable",
            df$resp,
            "Blocks",
            comp_x,
            comp_y,
            i_block,
            text = text,
            collapse =  collapse,
            no_Overlap = no_Overlap
        ) +
        geom_path(
            aes(x, y),
            data = circleFun(),
            col = "grey",
            size = 1
        ) +
        geom_path(
            aes(x, y),
            data = circleFun() / 2,
            col = "grey",
            size = 1,
            lty = 2
        )
    
    # remove legend if not on superblock
    if (!superblock || i_block != length(rgcca$a))
        p + theme(legend.position = "none")
    else
        p
}

#' Plot of components space
#'
#' Plots RGCCA components in a bi-dimensional space
#'
#' @inheritParams plotSamplesSpace
#' @inheritParams plotVariablesSpace
#' @param df A dataframe
#' @param title A character with the name of the space (either "Variables" or
#' "Samples")
#' @param group A vector of character with levels used to color the points
#' @param name_group A character giving the type of groups (either "Blocs" or
#' "Response")
#' @param p A ggplot object
#' @param colours A vectof of character to color quantitative dat
#' @examples
#' df = as.data.frame(matrix(runif(20*2, min = -1), 20, 2))
#' AVE = lapply(seq_len(4), function(x) runif(2))
#' rgcca.res = list(AVE = list(AVE_X = AVE))
#' plotSpace(rgcca.res, df, "Samples", rep(c("a","b"), each=10), "Response")
#' @export
plotSpace <- function(
    rgcca,
    df,
    title = "Biplot",
    group,
    name_group = "Response",
    comp_x = 1,
    comp_y = 2,
    i_block = NULL,
    p = NULL,
    text = TRUE,
    i_block_y = i_block,
    colours = c("blue", "gray", "#cd5b45"),
    collapse = FALSE,
    no_Overlap = FALSE) {

    if (!isTRUE(text)) {
        func <- quote(geom_point(size = PCH_TEXT_CEX))
        if (!is.numeric(na.omit(group)))
            func$mapping <- aes(shape = as.factor(group))
    } else {

        f <- "geom_text"
        func <- quote(
            get(f)(aes(label = rownames(df)),
            size = PCH_TEXT_CEX)
        )
        
        if (no_Overlap && nrow(df) <= 200) {
            f = paste0(f, '_repel') 
            func$force = 0.2
            func$max.iter = 500
        }
    }

    if (is.null(p))
        p <- ggplot(df, aes(df[, 1], df[, 2], colour = as.factor(group)))

    if (length(name_group) > 15)
        name_group <- name_group[seq_len(15)]

    if (is.null(name_group))
        name_group <- 0
    
    axis <- function(margin){
        element_text(
            face = "italic",
            size = AXIS_TITLE_CEX * 0.75,
            margin = margin
        )
    }

    p <- p + eval(as.call(func)) + theme_classic() + geom_vline(
            xintercept = 0,
            col = "grey",
            linetype = "dashed",
            size = 1
        ) + geom_hline(
            yintercept = 0,
            col = "grey",
            linetype = "dashed",
            size = 1
        ) + labs(
                title = paste(title, "space"),
                x = printAxis(rgcca, comp_x, i_block),
                y = printAxis(rgcca, comp_y, i_block_y),
            color = name_group,
            shape = name_group
        ) + 
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(breaks = NULL) +
        theme_perso() +
        theme(
            legend.key.width = unit(nchar(name_group), "mm"),
            axis.text = element_blank(),
            axis.title.y = axis(margin(0, 20, 0, 0)),
            axis.title.x = axis(margin(20, 0, 0, 0))
        )

    if (length(unique(group)) != 1 && title == "Variable") {
        orderColorPerBlocs(rgcca$a, p, collapse = collapse)
        # For qualitative response OR no response
    } else if ( isCharacter(group[!is.na(group)]) ||
                length(unique(group)) <= 5 || 
            all( levels(as.factor(group)) %in% c("obs", "pred") )
        ) {
        p + scale_color_manual(values = colorGroup(group))
        # quantitative response
    } else
        p + scale_color_gradientn(colours = colours, na.value = "black")

}

# Reorder the color of a ggplot (in case of a missing modality)
orderColorPerBlocs <- function(df, p, matched = NULL, collapse = FALSE) {

    J <- names(df)
    if (is.null(matched)) {
        matched <- seq_len(length(J))
        f <- "color"
    } else
        f <- "fill"

    if (collapse)
        n <- J
    else
        n <- J[-length(J)]

    func <- quote(get(paste("scale", f, "manual", sep = "_"))(
        values = colorGroup(J)[matched],
        limits = n[matched],
        labels = n[matched]
    ))
    
    return(p + eval(as.call(func)))
}

#' Histogram of a fingerprint
#'
#' Histogram of the higher outer weight vectors for a component of a block 
#' (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @inheritParams plotVariablesSpace
#' @param comp An integer giving the index of the analysis components
#' of a superblock
#' @param type A string giving the criterion to selects biomarkers : either 
#' "cor" for correlation between the component and the block
#' or "weight" for the weight of the RGCCA
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
# weights = lapply(seq_len(3), function(x) matrix(runif(7*2), 7, 2))
#' for(i in seq(3))
# row.names(weights[[i]]) <- paste0(letters[i],
#      letters[seq_len(nrow(weights[[i]]))])
# weights[[4]] = Reduce(rbind, weights)
# rgcca.res = list(a = weights)
# names(rgcca.res$a) = LETTERS[seq_len(4)]
# # With the 1rst component of the superblock
# plotFingerprint(rgcca.res, NULL, 1, TRUE, type = "weigth")
# # With the 2nd component of the 1rst block by selecting the ten higher weights
# plotFingerprint(rgcca.res, NULL, 2, FALSE, 10, 1, type = "weigth")
# library(RGCCA)
# data("Russett")
# blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#     politic = Russett[, 6:11] )
# rgcca.res = rgcca.analyze(blocks)
# getVariablesIndexes(rgcca.res, blocks, superblock = FALSE)
# df = getVariablesIndexes(rgcca.res, blocks, collapse = TRUE)
# plotFingerprint(rgcca.res, blocks, collapse = TRUE)
#' @export
plotFingerprint <- function(
    rgcca,
    blocks = NULL,
    comp = 1,
    superblock = TRUE,
    n_mark = 100,
    i_block = length(rgcca$a),
    type = "cor",
    collapse = FALSE) {

    df <- getVariablesIndexes(
        rgcca = rgcca,
        blocks = blocks,
        comp_x = comp,
        comp_y = comp,
        i_block = i_block,
        type = type,
        superblock = superblock,
        n_mark = n_mark,
        collapse = collapse,
        removeVariable = FALSE
    )
    
    J <- names(rgcca$a)

    title <- ifelse(type == "cor",
            "Variable correlations with",
            "Variable weights on")

    # sort in decreasing order
    df <- data.frame(getRankedValues(df, 1, TRUE), order = nrow(df):1)

    # max threshold for n
    if (nrow(df) >= n_mark)
        df <- df[seq(n_mark), ]

    # if the superblock is selected, color the text of the y-axis according
    # to their belonging to each blocks
    if (superblock & (collapse | (i_block == length(rgcca$a)))) {
        color <- factor(df$resp)
        levels(color) <- colorGroup(color)
        p <- ggplot(df, aes(order, df[, 1], fill = df$resp))
    } else {

        color <- "black"
        p <- ggplot(df, aes(order, df[, 1], fill = abs(df[, 1])))
    }

    p <- plotHistogram(p, df, title, as.character(color)) + 
        labs(subtitle = printAxis(rgcca, comp, i_block))

    # If some blocks have any variables in the top hit, selects the ones
    # corresponding
    if (collapse)
        col <- J
    else
        col <- J[-length(J)]

    matched <- match(rev(unique(df$resp)), col)

    # Force all the block names to appear on the legend
    if (length(color) != 1)
        p <- orderColorPerBlocs(rgcca$a, p, matched, collapse)
    if ( !superblock | (!collapse & i_block != length(rgcca$a)))
            p <- p + theme(legend.position = "none")

    return(p)
}

#' Histogram of Average Variance Explained
#'
#' Histogram of the model quality (based on Average Variance Explained)
#' for each blocks and sorted in decreasing order
#'
#' @inheritParams plotSamplesSpace
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' random_val = function(y=1) lapply(seq_len(4),
#' function(x) matrix(runif(4), y, 2))
#' rgcca.res = list(AVE = list(AVE_X = random_val()),
#'      a = random_val(2), ncomp = rep(2, 4))
#' names(rgcca.res$a) <- LETTERS[seq_len(4)]
#' library("ggplot2")
#' for(i in seq(1,4))
#' names(rgcca.res$AVE$AVE_X[[i]]) <- c(1,2)
#' plotAVE(rgcca.res)
#' @export
plotAVE <- function(rgcca) {

    ave <- 100 * unlist(rgcca$AVE$AVE_X)
    blocks <- factor(unlist(lapply(seq_len(length(names(rgcca$a))), 
            function(x) rep(names(rgcca$a)[x], rgcca$ncomp[x]))),
        levels = names(rgcca$a))
    ncomp <- as.factor(names(ave))

    y_ave_cum <- lapply(
        lapply(rgcca$AVE$AVE_X, 
            function(x) round(100 * cumsum(x), 1)), 
        function(x) c(0, x))
    y_ave_cum <- unlist(lapply(y_ave_cum, function(x)
            unlist(lapply(seq_len(length(x)),
                function(i) (x[i - 1] + x[i]) / 2))))

    ave_label <- unlist(lapply(rgcca$AVE$AVE_X, function(x)
            round(100 * x, 1)))
    ave_label[ave_label < max(y_ave_cum) / 20] <- ""

    df <- data.frame(ave, blocks, ncomp, stringsAsFactors = FALSE)

    p <- ggplot(data = df, 
        aes(
            x = blocks,
            y = ave,
            fill = ncomp,
            label = ave_label
        ))

    p <- plotHistogram(p, df, "Average Variance Explained") +
        scale_fill_manual(
            values = colorGroup(levels(df$ncomp)),
            labels = gsub("comp", " ", levels(df$ncomp))) + 
        geom_col(position = position_stack(reverse = TRUE)) +
        labs(subtitle = printAxis(rgcca, outer = TRUE)) +
        geom_text(aes(y = y_ave_cum),  cex = 3.5 * CEX, color = "white") +
        labs(fill = "Components")

    return(p)
}

#' Histogram settings
#'
#' Default font for a vertical barplot.
#'
#' @inheritParams plotSpace
#' @param color A vector of character giving the colors for the rows
#' @param low_col A character giving the color used for the lowest part of
#' the gradient
#' @param high_col A character giving the color used for the highest part of
#' the gradient
#' @param mid_col A character giving the color used for the middle part of
#' the gradient
#' @examples
#' df = data.frame(x = runif(30), order = 30:1)
#' library("ggplot2")
#' p = ggplot(df, aes(order, x))
#' plotHistogram(p, df, "This is my title")
#' # Add colors per levels of a variable
#' df$color = rep(c(1,2,3), each=10)
#' p = ggplot(df, aes(order, x, fill = color))
#' plotHistogram(p, df, "Histogram", as.character(df$color))
#' @export
plotHistogram <- function(
    p,
    df,
    title = "",
    color = "black",
    low_col = "khaki2",
    high_col = "coral3",
    mid_col = NULL) {
    

    if (nrow(df) <= 10 || title == "Average Variance Explained") {
        WIDTH <- NULL
        if (title == "Average Variance Explained")
            AXIS_TEXT_CEX <- 12
    } else
        WIDTH <- 1

    if (nrow(df) < 3)
        MAR <- 60
    else if (nrow(df) < 5)
        MAR <- 30
    else
        MAR <- 0

    axis <- function(margin){
        element_text(
            size = AXIS_TEXT_CEX,
            face = "italic",
            color = "gray40"
        )
    }

    p <- p + geom_bar(stat = "identity", width = WIDTH) +
        coord_flip() + labs(title = title,  x = "", y = "") +
        theme_classic() +
        theme_perso() +
        theme(
            axis.text.y = axis(),
            axis.text.x = axis(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.subtitle = element_text(
                hjust = 0.5,
                size = SUBTITLE_CEX,
                face = "italic"
            ),
            plot.margin = margin(0, 0, MAR, 0, "mm")
    )

    if (title != "Average Variance Explained") {
        p <- p +
            scale_x_continuous(breaks = df$order, labels = rownames(df)) +
            labs(fill = "Blocks")
        if (length(color) == 1) {
            if (is.null(mid_col))
                p <- p +
                    scale_fill_gradient(low = low_col, high = high_col) +
                    theme(legend.position = "none")
            else
                p <- p +
                    scale_fill_gradient2(low = low_col, high = high_col, mid = mid_col)
        }
    }

    return(p)
}

#' Variable contribution
#' 
#' Extract the contibution of variables to the model by using correlation or weight
#' @inheritParams plotVariablesSpace
#' @param type A character giving the choice ot the index between cor or weight
#' @param comp_z An integer giving the index of the analysis component used
#' for the z-axis
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq_len(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks, ncomp = c(3,2,4))
#' getVar(rgcca.res, blocks)
#' # On the first block and with weights
#' getVar(rgcca.res, blocks, 2, 1, i_block = 1, type = "weights")
#' # With 3 components and on the variables of two blocks
#' superblocks <- rep(list(Reduce(cbind, c(blocks[1], blocks[3]))), 2)
#' names(superblocks) <- names(blocks)[c(1, 3)]
#' rgcca.res = rgcca.analyze(blocks[c(1,3)], ncomp = c(3,4))
#' getVar(rgcca.res, superblocks, comp_z = 3, i_block = 1, type = "cor", collapse = TRUE)
#' getVar(rgcca.res, superblocks, 2, 1, 3, 1, "weights", TRUE)
#' @return A dataframe containing the indexes for each selected components
#' @export
getVar <- function(
    rgcca,
    blocks = NULL,
    comp_x = 1,
    comp_y = 2,
    comp_z = NULL,
    i_block = length(blocks),
    type = "cor",
    collapse = FALSE) {

    if (!collapse)
        i_block_2 <- i_block
    else
        i_block_2 <- 1

    if (is.null(blocks))
        row.names = row.names(rgcca$a[[i_block]])
    else
        row.names = colnames(blocks[[i_block]])

    if (type == "cor")
        f <- function(x) cor(
                blocks[[i_block_2]],
                rgcca$Y[[i_block]][, x],
                use = "pairwise.complete.obs"
            )
    else{
        if (!collapse)
            f <- function(x) rgcca$a[[i_block]][, x]
        else
            f <- function(x) unlist(
                sapply(
                    1:length(blocks),
                    function(y) rgcca$a[[y]][, x]
                )
            )
    }

    data.frame(
        sapply(
            c(comp_x, comp_y, comp_z[comp_z >= rgcca$ncomp[i_block]]), 
            function(x) f(x)
        ),
        row.names = row.names
    )

}

#' Rank values of a dataframe in decreasing order
#'
#' @param df A dataframe
#' @param comp An integer giving the index of the analysis components
#' @param allCol A boolean to use all the column of the dataframe
#' @return A datafram with ordered values
#' @examples 
#' df = sapply(seq(2), function(x) runif(10))
#' getRankedValues(df)
#' @export
getRankedValues <- function(df, comp = 1, allCol = TRUE) {
    
    ordered <- order(abs(df[, comp]), decreasing = TRUE)

    if (allCol)
        comp <- seq_len(ncol(df))

    res <- df[ordered, comp]

    if (!allCol)
        names(res) <- row.names(df)[ordered]

    return(res)
}


# Print and save variables analysis attributes
saveVars <- function(
    rgcca,
    blocks,
    comp_x = 1,
    comp_y = 2,
    file = "variables.tsv") {

    indexes <- c("cor", "weight")

    vars <- Reduce(rbind, lapply(seq_len(length(blocks)), function(i)
            data.frame(
                Reduce(cbind,
                        lapply(indexes, function(x)
                            getVar(rgcca, blocks, comp_x, comp_y, i_block = i, type = x))),
                names(blocks)[i]
            )))

    colnames(vars) <- c(as.vector(sapply(indexes, function(x)
            paste0(x, ".", paste0("axis.", c(comp_x, comp_y))))), "block")

    write.table(vars, file, sep = "\t")

    invisible(vars)
}

# Print and save indidvidual analysis attributes
saveInds <- function(
    rgcca,
    blocks,
    comp_x = 1,
    comp_y = 2,
    file = "individuals.tsv") {
    
    inds <- Reduce(cbind, lapply(
        rgcca$Y,
        function(x) x[, c(comp_x, comp_y)]))
    colnames(inds) <- as.vector(sapply(
        names(blocks),
        function(x) paste0(x, ".axis", c(comp_x, comp_y))))

    write.table(inds, file, sep = "\t")

    invisible(inds)
}
