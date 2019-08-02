# Author: Etienne CAMENEN
# Date: 2018
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: A user-friendly multi-blocks analysis (Regularized Generalized Canonical Correlation Analysis, RGCCA)
# with all default settings predefined. Produce four figures to help clinicians to identify fingerprint:
# the samples and the variables projected on the two first component of the multi-block analysis, the histograms
# of the most explicative variables and the explained variance for each blocks.

#Global settings
AXIS_TITLE_SIZE = 19
AXIS_TEXT_SIZE = 10
PCH_TEXT_SIZE = 3
AXIS_FONT = "italic"
SAMPLES_COL_DEFAULT = "brown3"

# X- Y axis format for plotly objets : no axis, no ticks
ax <- list(linecolor = "white", ticks = "", titlefont = list(size = 23))
ax2 <- list(linecolor = "white", tickfont = list(size = 10, color = "grey"))

# Dynamic visualization of the outputs
# f: ggplot2 function
# ax: list object containing attributes of xaxis / yaxis parameter in plotly (https://plot.ly/javascript/reference/, xaxis/yaxis)
# text: axis information to print (among y, x, x+y, text) (https://plot.ly/javascript/reference/, hoverinfo)
# dynamicTicks: a bolean giving the generation of axis tick labels (otherwhise samplesPlot which do not have traces could not be convereted in ggplotly)
# return a plotly object
dynamicPlot = function (f, ax, text = "name+x+y", legend = TRUE, dynamicTicks = FALSE) {

  # Convert a ggplot into a plotly object
  # add a layout with predefined formats for x- and y- axis
  # set the style to show onMouseOver text
  p = plotly_build( ggplotly(f, dynamicTicks = dynamicTicks) %>%
                     layout(xaxis = ax, yaxis = ax, annotations = list(showarrow = F, text = "")) %>%
                     style(hoverinfo = text))

  if(legend){
    # set on the top the position of the legend title
    p$x$layout$annotations[[1]]$yanchor = "top"

    # Deals with a too short name of modalities
    p$x$layout$margin$r = nchar(p$x$layout$annotations[[1]]$text) * 13
    p$x$layout$margin$t = 100

    # for shiny corcircle, if text = TRUE, two legends will appear. only the first one will be selected
    title = unlist(strsplit(p$x$layout$annotations[[1]]$text, "<br />"))[1]
    # to prevent print a "NA" when there is no legend in plot
    if(is.na(title))
      title = ""
    # set the font for this title
    p$x$layout$annotations[[1]]$text = paste0("<i><b>", title, "</b></i>")

    #Sys.info()[['sysname']]

    if(!is.null(f$labels$subtitle)){
       if(packageVersion("plotly") < 4.9)
         p$x$layout$title = paste0(p$x$layout$title, '<br><b>', "c" , substring(f$labels$subtitle, 2), '</b>')
      else
        p$x$layout$title$text = paste0(p$x$layout$title$text, '<br><b>', "c" , substring(f$labels$subtitle, 2), '</b>')
    }
  }

  if( ncol(f$data) == 3 )
    p$sample_names = lapply(levels(as.factor(f$data[, 3])), function(x) row.names(subset(f$data, f$data[, 3] == x)))
  else
    p$sample_names = list(row.names(f$data))

  return(p)
}


dynamicPlotBoot = function(p){

  p = dynamicPlot(p, ax2, "text")
  n = length(p$x$data)
  m = unlist(lapply(p$x$data, function(x) !is.null(x$orientation)))
  j = length(m[m])

  for (i in 1:j){
    p$x$data[[i]]$text = paste( round(p$x$data[[i]]$x, 3), "+/-", round(p$x$data[[n]]$error_x$array[j], 3) )
    j = j - 1
  }

  # Remove the onMouseOver for the last trace
  changeText ( p ) %>%
    style(error_x = list( array = p$x$data[[n]]$error_x$array, color = "gray"), hoverinfo = "none", traces = n)
}

# p: a ggplot function
# hovertext : a boolean for the use of hovertext (if TRUE) as the attribute to parse the onMouseOver text ("text" attribute, if FALSE)
changeHovertext = function(p, hovertext = TRUE){

  attr = ifelse(hovertext, "hovertext", "text")
  # identify the order / id of the traces which corresponds to x- and y-axis (should be before the splitting function)
  traces = which(lapply(p$x$data, function(x) length(grep("intercept", x$text)) == 1) == T)
  # length of groups of points without traces and circle points
  n = which(sapply(p$x$data, function(x) match("xintercept: 0", x$text)==1)) -1

  for (i in 1:n) {

      # For each lines of each group in the legend
      for (j in 1:length(p$x$data[[i]][attr][[1]])){
        # Distinguish each duplicate by splitting with "<br>"  and separe them in key/value by splitting with ": " (like a dictionnary)
        l_text  = lapply( as.list( strsplit( p$x$data[[i]][attr][[1]][j], "<br />" )[[1]] ), function(x) strsplit(x, ": ")[[1]] )
        # keep only the (x, y) coordinates with the key df[, ] and the response if exists
        l_text = unlist(lapply(l_text, function(x, y) {
          if(x[1] %in% paste0("df[, ", c(1, 2), "]"))
            round(as.numeric(x[2]), 3)
          else if(x[1] == "resp")
            x[2]
          })
        )

        # do not print names because text = FALSE plots have'nt names
        name = ifelse(hovertext, paste0("name: ", p$x$data[[i]]$text[j], "<br />"), "")
        # Overwrite the onMouseOver text with the (x, y) coordinates and the response if exists
        p$x$data[[i]][attr][[1]][j] = paste0(name, "x: ", l_text[1], "<br />y: ", l_text[2], ifelse(length(l_text)==3,  paste0("<br />response: ", l_text[3]) , ""))

      }
  }
  # Remove the x- and y- axis onOverMouse
  ( style(p, hoverinfo = "none", traces = traces ) )
}

changeText = function(p){
  for (i in 1:length(p$x$data))
    p$x$data[[i]]$text = sub( "order: .*<br />df\\[, 1\\]: (.*)<.*", "\\1\\", p$x$data[[i]]$text)

  return (p)
}

# Creates a circle
circleFun = function(center = c(0, 0), diameter = 2, npoints = 100) {

  r = diameter/2
  tt = seq(0, 2 * pi, length.out = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[2] + r * sin(tt)

  return(data.frame(x = xx, y = yy))
}

#' Print the variance of a component
#'
#' Prints the percent of explained variance for a component of a block (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param n An integer giving the index of the analysis component
#' @param i An integer giving the index of a list of blocks
#' @param outer A boolean for ave plot case
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
#' @export printAxis
printAxis = function (rgcca, n = NULL, i = NULL, outer = FALSE){

  # by default, take the last block
  if ( is.null(i) )
    i = length(rgcca$AVE$AVE_X)

  nvar = varSelected(rgcca, i, n)

  if(class(rgcca) != "sgcca" | nvar == length(rgcca$a[[i]][, n]) )
    varText = ""
  else
    varText = paste0(nvar, " variables, ")

  ave = quote(paste0(round(AVE[n] * 100 , 1),"%"))

  if(isTRUE(outer)){
    AVE = rgcca$AVE$AVE_outer
    n = c(1, 2)
    paste0("First outer comp. : ", paste(eval(ave), collapse=" & "))
  }else{
    AVE = rgcca$AVE$AVE_X[[i]]
    paste0("Component ", n, " (", varText, eval(ave), ")")
  }
}

varSelected = function(rgcca, i_block, comp){
  # Get the variables with a weight != 0
  sum(rgcca$a[[i_block]][,comp] != 0)
}

#' Default font for plots
theme_perso = function() {

  theme(
    legend.text = element_text(size = 13),
    legend.title = element_text(face="bold.italic", size=16),
    plot.title = element_text(size = 25, face = "bold", hjust=0.5, margin = margin(0,0,20,0))
  )
}

colorGroup = function(group){
  palette = rep(c("#cd5b45", "#71ad65", "#3c78b4", "#ffc600",
                  "#b448af", "#9d9d9d", "#abcaef", "#4a6f43",  "#f0e500",
                  "#efb8f0", "black", "#d6d6d6" ), 10)
  palette[0 : length(levels(as.factor(group))) ]
}

#' Plot of samples space
#'
#' Plots samples on two components of a RGCCA
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param resp A vector of characters corresponding either to a qualitative variable with levels or a continuous variable
#' @param comp_x An integer giving the index of the analysis component used for the x-axis
#' @param comp_y An integer giving the index of the analysis component used for the y-axis
#' @param i_block An integer giving the index of a list of blocks
#' @param text A bolean to represent the points with their row names (TRUE) or with circles (FALSE)
#' @param i_block_y An integer giving the index of a list of blocks (another one, different from the one used in i_block)
#' @param reponse_name A character giving the legend title
#' @examples
#' coord = lapply(1:3, function(x) matrix(runif(15 * 2, min = -1), 15, 2))
#' AVE_X = lapply(1:3, function(x) runif(2))
#' rgcca.res = list(Y = coord, AVE = list(AVE_X = AVE_X))
#' # Using a superblock
#' plotSamplesSpace(rgcca.res, rep(LETTERS[1:3], each = 5))
#' # Using the first block
#' plotSamplesSpace(rgcca.res, runif(15, min=-15, max = 15), 1, 2, 1)
#' @export plotSamplesSpace
plotSamplesSpace = function (rgcca, resp, comp_x = 1, comp_y = 2, i_block = NULL, text = TRUE, i_block_y = NULL, reponse_name = "Response"){
  # resp : color the points with a vector

  # Avoid random with ggrepel
  set.seed(1)

  if ( is.null(i_block) )
    i_block = length(rgcca$Y)

  if(is.null(i_block_y))
    df = data.frame(rgcca$Y[[i_block]][, c(comp_x, comp_y)])
  else
    df =  data.frame(rgcca$Y[[i_block]][, comp_x], rgcca$Y[[i_block_y]][, comp_y] )

  if(nrow(df) > 100)
    PCH_TEXT_SIZE = 2

  # if the resp is numeric
  if (  length(unique(as.matrix(resp))) > 1 ){

    names = row.names(resp)
    resp = apply(as.matrix(resp), 1, as.character)

    if(!is.null(names)){

      resp = as.matrix(resp, row.names = names)
      name_blocks = row.names(blocks[[i_block]])
      diff_column = setdiff(name_blocks, names)

      if ( identical(diff_column, name_blocks ) ){
        warning("No match has been found with the row names of the group file.")
        resp <- rep("NA", nrow(df))
      }else{
        if(length(diff_column) > 0 ){
          resp[diff_column ] <- "NA"
          names(resp)[names(resp) == ""] <- names
        }else{
          names(resp) = names
        }
      resp = resp[row.names(blocks[[i_block]])]

      }
    }else{
      warning("No row names have been found in the group file.")
      resp <- rep("NA", nrow(df))
    }

    if( ! unique(isCharacter(as.vector(resp))) && length(unique(resp)) > 5 ){

      resp[resp == "NA"] <- NA
      resp = as.numeric(resp)
      df$resp = resp

      # add some transparency
      p = ggplot(df, aes(df[, 1], df[, 2], color =  resp))

    }else
      p = NULL
  }else
    p = NULL

  p = plotSpace(rgcca, df, "Sample", resp, reponse_name, comp_x, comp_y, i_block, p, text, i_block_y)

  # remove legend if missing
  if ( length(unique(resp)) == 1){
    p + theme(legend.position = "none")
  }else
    p
}

#' Get the blocs of each variables
#'
#' Get a vector of block names for each corresponding variable. The last block is considered as the superblock and ignored.
#'
#' @param df A list of matrix where their names are those of the blocks and the superblock and their rows are named after their variables
#' @return A vector of character giving block names for each corresponding variable.
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' rgcca.res = list(a = rep(NA, 4))
#' names(rgcca.res$a) = LETTERS[1:4]
#' getBlocsVariables(rgcca.res)
#' # a, b, c
#' @export getBlocsVariables
getBlocsVariables = function(df){

  rep( names(df)[-length(df)],
       sapply(df[1:(length(df)-1)], NROW))
}

#' Plot of variables space
#'
#' Correlation circle highlighting the contribution of each variables to the construction of the RGCCA components
#' @param rgcca A list giving the results of a R/SGCCA
#' @param blocks A list of matrix
#' @param comp_x An integer giving the index of the analysis component used for the x-axis
#' @param comp_y An integer giving the index of the analysis component used for the y-axis
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of a superblock
#' @param i_block An integer giving the index of a list of blocks
#' @param text A bolean to represent the points with their row names (TRUE) or with circles
#' @param removeVariable A bolean to keep only the 100 variables of each component with the biggest correlation#'
#' @param n_mark An integer giving the number of top variables to select
#' @examples
#' setMatrix = function(nrow, ncol, iter = 3) lapply(1:iter, function(x) matrix(runif(nrow * ncol), nrow, ncol))
#' blocks = setMatrix(10, 5)
#' blocks[[4]] = Reduce(cbind, blocks)
#' for (i in 1:4)
#'     colnames(blocks[[i]]) = paste0( LETTERS[i], as.character(1:NCOL(blocks[[i]])))
#' coord = setMatrix(10, 2, 4)
#' a = setMatrix(5, 2)
#' a[[4]] = matrix(runif(15 * 2), 15, 2)
#' AVE_X = lapply(1:4, function(x) runif(2))
#' rgcca.res = list(Y = coord, a = a, AVE = list(AVE_X = AVE_X))
#' names(rgcca.res$a) = LETTERS[1:4]
#' # Using a superblock
#' plotVariablesSpace(rgcca.res, blocks, 1, 2, TRUE)
#' # Using the first block
#' plotVariablesSpace(rgcca.res, blocks, 1, 2, FALSE, 1)
#' @export plotVariablesSpace
plotVariablesSpace = function(rgcca, blocks, comp_x = 1, comp_y = 2, superblock = TRUE, i_block = NULL, text = TRUE,
                              removeVariable = TRUE, n_mark = 100){

  x = y = selectedVar = NULL

  if ( is.null(i_block) )
    i_block = length(blocks)

  df =  getVar(rgcca, blocks, comp_x, comp_y, i_block, "cor")
  df_temp <- df

  if(class(rgcca)=="sgcca"){
    selectedVar = rgcca$a[[i_block]][,comp_x] != 0 | rgcca$a[[i_block]][,comp_y] != 0
    df = df[selectedVar, ]
  }else{

    if(n_mark < 2)
      n_mark = nrow(df)

    if(removeVariable & nrow(df) > 2 * n_mark){
      selectedVar = unique( as.vector (unique( sapply(c(1, 2),
                                         function(x) row.names(data.frame(df[order(abs(df[, x]), decreasing = TRUE),])[1:n_mark,])))))
  	  df = df[selectedVar, ]
    }
  }

  # if superblock is selected, color by blocks
  if ( superblock & ( i_block == length(blocks)) ){

    color = getBlocsVariables(rgcca$a)

    if(class(rgcca)=="sgcca"){
      names(color) = row.names(df_temp)
      color = color[row.names(df)]
    }else{
    	if(!is.null(selectedVar))
    		color = color[unlist(lapply(1:length(selectedVar), function(x) which(colnames(blocks[[length(blocks)]]) == selectedVar[x])))]
    }

  }else{
    color = rep(1, NROW(df))
  }

  df = data.frame(df, color)


  p = plotSpace(rgcca, df, "Variable", color, "Blocks", comp_x, comp_y, i_block, text = text) +
    geom_path(aes(x, y), data = circleFun(), col = "grey", size = 1) +
    geom_path(aes(x, y), data = circleFun()/2, col = "grey", size = 1, lty = 2)

  # remove legend if not on superblock
  if ( !superblock || i_block != length(blocks) )
    p + theme(legend.position = "none")
  else
    p

}

#' Plot of components space
#'
#' Plots RGCCA components in a bi-dimensional space
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param df A dataframe
#' @param title A character with the name of the space (either "Variables" or "Samples")
#' @param group A vector of character with levels used to color the points
#' @param name_group A character giving the type of groups (either "Blocs"  or "Response")
#' @param comp_x An integer giving the index of the analysis component used for the x-axis
#' @param comp_y An integer giving the index of the analysis component used for the y-axis
#' @param i_block An integer giving the index of a list of blocks
#' @param p A ggplot object
#' @param text A bolean to represent the points with their row names (TRUE) or with circles (FALSE)
#' @param i_block_y An integer giving the index of a list of blocks (another one, different from the one used in i_block)
#' @param colours A vectof of character to color quantitative data
#' @examples
#' df = as.data.frame(matrix(runif(20*2, min = -1), 20, 2))
#' AVE =  lapply(1:4, function(x) runif(2))
#' rgcca.res = list(AVE = list(AVE_X = AVE))
#' plotSpace(rgcca.res, df, "Samples", rep(c("a","b"), each=10), "Response")
#' @export plotSpace
plotSpace = function (rgcca, df, title, group, name_group, comp_x = 1, comp_y = 2, i_block = NULL, p = NULL, text = TRUE,
                      i_block_y = NULL, colours = c("blue", "gray", SAMPLES_COL_DEFAULT)){

  if(is.null(i_block_y))
    i_block_y = i_block

  if (!isTRUE(text)){
    func = quote(geom_point(size = PCH_TEXT_SIZE))
    if (!is.numeric(na.omit(group)))
      func$mapping = aes(shape = as.factor(group))
  }else{
    f = "geom_text"
    func = quote(get(f)(aes(label = rownames(df)), size = PCH_TEXT_SIZE))

    # if(no_Overlap && nrow(df) <= 100){
    #   f = paste0(f, "_repel")
    #   func$force = 0.2
    #   func$max.iter = 500
    # }
  }

  if (is.null(p)){
    p = ggplot(df, aes(df[,1], df[,2], colour = as.factor(group)))
  }

  if(length(name_group) > 15)
    name_group <- name_group[1:15]

  if(is.null(name_group))
    name_group <- 0

  p = p + eval(as.call(func)) +
    theme_classic() +
    geom_vline(xintercept = 0, col = "grey", linetype = "dashed", size = 1) +
    geom_hline(yintercept = 0, col = "grey", linetype = "dashed", size = 1) +
    labs ( title = paste(title, "space"),
  		 x = printAxis(rgcca, comp_x, i_block),
  		 y = printAxis(rgcca, comp_y, i_block_y),
       color = name_group,
       shape = name_group) +
    scale_y_continuous(breaks = NULL) +
    scale_x_continuous(breaks = NULL) +
    theme_perso() +
    theme(
      legend.key.width = unit(nchar(name_group), "mm"),
      axis.text = element_blank(),
      axis.title.y = element_text(face = AXIS_FONT, margin = margin(0,20,0,0), size = AXIS_TITLE_SIZE),
      axis.title.x = element_text(face = AXIS_FONT, margin = margin(20,0,0,0), size = AXIS_TITLE_SIZE)
    )

  if(length(unique(group)) != 1 && title == "Variable"){
    orderColorPerBlocs(rgcca, p)
  # For qualitative response OR no response
  }else if(isCharacter(group[!is.na(group)]) || length(unique(group)) <= 5 ){
    p + scale_color_manual(values = colorGroup(group))
  # For quantitative response
  }else
    p + scale_color_gradientn(colours = colours, na.value = "black")

}

orderColorPerBlocs <- function(rgcca, p, matched  = NULL){

  J <- names(rgcca$a)

  if(is.null(matched)){
    matched <- 1:length(J)
    f <- "color"
  }else
    f <- "fill"

  func <- quote(get(paste("scale", f, "manual", sep = "_"))(values = colorGroup(J)[matched],
                                                           limits = J[-length(J)][matched],
                                                           labels = J[-length(J)][matched]))

  return( p + eval(as.call(func)) )
}

#' Histogram of a fingerprint
#'
#' Histogram of the higher outer weight vectors for a component of a block (by default, the superblock or the last one) analysed by R/SGCCA
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @param blocks A list of matrix
#' @param comp An integer giving the index of the analysis components
#' @param n_mark An integer giving the number of top variables to select
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of a superblock
#' @param i_block An integer giving the index of a list of blocks
#' @param type A string giving the criterion to selects biomarkers : either "cor" for correlation between the component and the block
#' or "weight" for the weight of the RGCCA
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' weights = lapply(1:3, function(x) matrix(runif(7*2), 7, 2))
#' weights[[4]] = Reduce(rbind, weights)
#' rgcca.res = list(a = weights)
#' names(rgcca.res$a) = LETTERS[1:4]
#' for(i in seq(1,4))
#' row.names(rgcca.res$a[[i]]) <- paste0(letters[i],letters[1:nrow(rgcca.res$a[[i]])])
#' # With the 1rst component of the superblock
#' plotFingerprint(rgcca.res, NULL, 1, TRUE, type = "weigth")
#' # With the 2nd component of the 1rst block by selecting the ten higher weights
#' plotFingerprint(rgcca.res, NULL, 2, FALSE, 10, 1, type = "weigth")
#' @export plotFingerprint
plotFingerprint = function(rgcca, blocks = NULL, comp = 1, superblock = TRUE, n_mark = 100, i_block = NULL, type = "cor"){

  color = NULL
  J = names(rgcca$a)

  # if no specific block is selected, by default, the superblock is selected (or the last one)
  if ( is.null(i_block) )
    i_block = length(rgcca$a)

  title = ifelse(type == "cor", "Variable correlations with", "Variable weights on")
  criterion = getVar(rgcca, blocks, comp, comp, i_block, type)

  # select the weights (var to add a column to work with comp = 1)
  df = data.frame(criterion, var = row.names(criterion))

  # Get a qualitative variable with which block is associated with each variables
  if (  superblock & ( i_block == length(rgcca$a) ) )
    df = data.frame( df, color = as.factor(getBlocsVariables(rgcca$a)) )

  # sort in decreasing order
  df = data.frame(getRankedValues(df, 1, T), order = nrow(df):1)

  # selected variables in sgcca
  nvar_select = varSelected(rgcca, i_block, comp)
  if(n_mark > nvar_select)
    n_mark = nvar_select

  # max threshold for n
  if(NROW(df) >= n_mark)
    df = df[1:n_mark,]

  # if the superblock is selected, color the text of the y-axis according to their belonging to each blocks
  if (  superblock & ( i_block == length(rgcca$a) ) ){
    color2 = factor(df$color); levels(color2) = colorGroup(color2)
  }else{
    color2 = "black"
  }

  if (  superblock & i_block == length(rgcca$a) ){
    # levels(df$color) = rev(levels(df$color))
    p = ggplot(df, aes(order, df[, 1], fill = color))
  }else{
    p = ggplot(df, aes(order, df[, 1], fill = abs(df[, 1])))
  }

  p = plotHistogram(p, df, title, as.character(color2)) +
    labs(subtitle = printAxis(rgcca, comp, i_block))

  # If some blocks have any variables in the top hit, selects the ones corresponding
  matched <- match(rev(unique(df$color)), J[-length(J)])

  # Force all the block names to appear on the legend
  if(length(color2) != 1)
    p = orderColorPerBlocs(rgcca, p, matched)
  if (  !superblock | i_block != length(rgcca$a) )
    p = p + theme(legend.position = "none")

  return(p)
}

#' Histogram of Average Variance Explained
#'
#' Histogram of the model quality (based on Average Variance Explained) for each blocks and sorted in decreasing order
#'
#' @param rgcca A list giving the results of a R/SGCCA
#' @seealso \code{\link[RGCCA]{rgcca}}, \code{\link[RGCCA]{sgcca}}
#' @examples
#' random_val = function(y=1) lapply(1:4, function(x) matrix(runif(4), y, 2))
#' rgcca.res = list(AVE = list(AVE_X = random_val()), a = random_val(2), ncomp = rep(2, 4))
#' names(rgcca.res$a) <- LETTERS[1:4]
#' library("ggplot2")
#' for(i in seq(1,4))
#' names(rgcca.res$AVE$AVE_X[[i]]) <- c(1,2)
#' plotAVE(rgcca.res)
#' @export plotAVE
plotAVE = function(rgcca){

  ave = 100 * unlist(rgcca$AVE$AVE_X)
  blocks = factor(unlist(lapply(1:length(names(rgcca$a)), function(x) rep(names(rgcca$a)[x], rgcca$ncomp[x]))), levels = names(rgcca$a))
  ncomp = as.factor(names(ave))

  y_ave_cum = lapply(lapply(rgcca$AVE$AVE_X, function(x) round(100 * cumsum(x), 1)), function(x) c(0, x))
  y_ave_cum = unlist(lapply(y_ave_cum, function(x) unlist(lapply(1:length(x), function(i) (x[i-1] + x[i]) / 2 ))))

  ave_label = unlist(lapply(rgcca$AVE$AVE_X, function(x) round(100 * x, 1)))
  ave_label[ave_label < max(y_ave_cum)/20] =  ""

  df = data.frame(ave, blocks, ncomp, stringsAsFactors = F)
  p = ggplot(data=df, aes(x=blocks, y=ave, fill = ncomp, label =  ave_label))

  p = plotHistogram(p, df, "Average Variance Explained") +
    scale_fill_manual(values=colorGroup(levels(df$ncomp)), labels = gsub("comp", " ", levels(df$ncomp))) +
    geom_col(position = position_stack(reverse = TRUE)) +
    labs(subtitle = printAxis(rgcca, outer = TRUE)) +
    geom_text(aes(y = y_ave_cum), cex = 3.5, color = "white") +
    labs( fill = "Components" )

  return(p)
}

#' Histogram settings
#'
#' Default font for a vertical barplot.
#'
#' @param p A ggplot object.
#' @param df A dataframe with a column named "order"
#' @param title A character string giving a graphic title
#' @param color A vector of character giving the colors for the rows
#' @param low_col A character giving the color used for the lowest part of the gradient
#' @param high_col A character giving the color used for the highest part of the gradient
#'
#' @examples
#' df = data.frame(x = runif(30), order = 30:1)
#' library("ggplot2")
#' p = ggplot(df, aes(order, x))
#' plotHistogram(p, df, "This is my title")
#' # Add colors per levels of a variable
#' df$color = rep(c(1,2,3), each=10)
#' p = ggplot(df, aes(order, x, fill = color))
#' plotHistogram(p, df, "Histogram", as.character(df$color))
#' @export plotHistogram
plotHistogram = function(p, df, title = "", color = "black", low_col = "khaki2", high_col = "coral3"){

  if ( nrow(df) <= 10 || title == "Average Variance Explained" ){
    WIDTH = NULL
    if(title == "Average Variance Explained")
      AXIS_TEXT_SIZE = 12
  }else
    WIDTH = 1

  if(nrow(df) < 3)
    MAR = 60
  else if(nrow(df) < 5)
    MAR = 30
  else
    MAR = 0

  p = p +
    geom_bar(stat = "identity", width = WIDTH) +
    coord_flip()  +
    labs(
      title = title,
      x = "", y = "") +
    theme_classic() +
    theme_perso() +
    theme(
      axis.text.y = element_text(size = AXIS_TEXT_SIZE, face = AXIS_FONT, color = "gray40"),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE, face = AXIS_FONT, color = "gray40"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"),
      plot.margin = margin(0, 0, MAR, 0, "mm"))

  if(title != "Average Variance Explained"){
    p  = p +
      scale_x_continuous(breaks = df$order, labels = rownames(df)) +
      labs( fill = "Blocks")
    if(length(color) == 1){
      p = p + scale_fill_gradient(low = low_col, high = high_col) +
      theme(legend.position = "none")
    }
  }

  return(p)
}

 corResponse = function(rgcca, blocks, i_response = NULL, comp = 1, i_block = 1){

  if(is.null(i_response))
    response = blocks[[ length(rgcca$a) ]]
  else{
    response = blocks[[i_response]]
    diff_column = setdiff(row.names(blocks[[i_block]]), row.names(response))
    response[diff_column, ] <- NA
    response = response[row.names(rgcca$Y[[i_block]]),]
    #options(warn = -1)
    # Disabling automatic factor conversion for some columns
    response = apply(response, 2, as.double)
    #options(warn = 0)
  }

  cor.res = matrix(cor(rgcca$Y[[i_block]][, comp],
                       response,
                       use = "pairwise.complete.obs"),
                   dimnames = list(colnames(response, NULL)))
  res = data.frame(cor = cor.res[order(abs(cor.res[, 1]),
                                       decreasing = TRUE),],
                   order = length(cor.res):1)

  if(VERBOSE) print(res)

  p = ggplot(res, aes(order, cor, fill = abs(res[,1])))
  plotHistogram(p, res, "Correlation with response") +
    labs(subtitle = printAxis(rgcca, comp, i_block))  +
    theme(legend.position = "none")
}

# blocks : should not be null with cor mode
getVar = function(rgcca, blocks = NULL, comp_x = 1, comp_y = 2, i_block = NULL, type = "cor"){

  if ( is.null(i_block) )
    i_block = length(blocks)

  if(type == "cor")
    f = function(x) cor( blocks[[i_block]], rgcca$Y[[i_block]][, x] )
  else
    f = function(x) rgcca$a[[i_block]][, x]

  return(  data.frame(
    #correlation matrix within a block for each variables and each component selected
    sapply ( c(comp_x, comp_y), function(x) f(x)) ,
    row.names = row.names(rgcca$a[[i_block]])
  ))

}

#' Rank values of a dataframe in decreasing order
#'
#' @param df A dataframe
#' @param comp An integer giving the index of the analysis components
#' @param allCol A boolean to use all the column of the dataframe
getRankedValues = function(df, comp = 1, allCol = T){

  ordered = order(abs(df[, comp]), decreasing = T)

  if (allCol)
    comp = 1:ncol(df)

  res = df[ordered, comp]

  if (!allCol)
    names(res) = row.names(df)[ordered]

  return(res)
}


# Print variables analysis attributes
saveVars <- function(rgcca, blocks, comp_x = 1, comp_y= 2, file = "variables.tsv"){

  indexes <- c("cor", "weight")

  vars <- Reduce(rbind, lapply( 1: length(blocks), function(i)
    data.frame(Reduce(cbind, lapply(indexes, function(x) getVar(rgcca, blocks, comp_x, comp_y, i, x) )), names(blocks)[i])
  ))

  colnames(vars) <- c(as.vector(sapply(indexes, function(x) paste0(x, ".", paste0("axis.", c(comp_x, comp_y))))), "block")
  write.table(vars, file, sep = "\t")
  invisible(vars)
}

saveInds <- function(rgcca, blocks, comp_x = 1, comp_y = 2, file = "individuals.tsv"){

  inds <- Reduce(cbind,lapply(rgcca$Y, function(x) x[, c(comp_x, comp_y)]))
  colnames(inds) <- as.vector(sapply(names(blocks), function(x) paste0(x, ".axis", c(comp_x, comp_y))))
  write.table(inds, file, sep = "\t")
  invisible(inds)
}
