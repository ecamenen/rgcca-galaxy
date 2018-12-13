circleFun = function(center = c(0, 0), diameter = 2, npoints = 100) {
  # Creates x, y coordinates for a circle

  r = diameter/2
  tt = seq(0, 2 * pi, length.out = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

printAxis = function (rgcca, n, i = NULL){
  # Prints the % of explained variance for an axis
  # n: number of the axis
  # i: index of the blocks

  if ( is.null(i) )
    i = length(blocks)
  paste("Axis ", n, " (", round(rgcca$AVE$AVE_X[[i]][n] * 100 , 1),"%)", sep="")
}

theme_perso = function() {
  # Default font for plots

  theme(
    legend.text = element_text(size = 13),
    legend.title = element_text(face="bold.italic", size=16),
    plot.title = element_text(size = 25, face = "bold", hjust=0.5, margin = margin(0,0,20,0))
  )
}

plotSamplesSpace = function (rgcca, compX, compY, i_block=NULL, group=NULL){
  # Projectes coordinates of samples in a bi-dimensional space
  # compX: component used for the x-axis
  # compY: component used for the y-axis
  # i_block : index of the block
  # group : color the points with a vector

  if ( is.null(i_block) )
    i_block = length(blocks)

  df = data.frame(rgcca$Y[[i_block]])

  # if the group is numeric
  if ( ! is.null(group) ){
    if( ! unique(isCharacter(as.vector(group)))){
      # add some transparency
      p = ggplot(df, aes(df[,compX], df[,compY], alpha = (group - min(group)) / max(group - min(group)))) +
        # get a color scale by quantile
            scale_alpha_continuous(
            name = "group",
            breaks = seq(0,1,.25),
            labels = round(quantile(group))
        ) + geom_text(color = COLOR_SAMPLES_DEF, aes(label = rownames(df)), size = PCH_TEXT_SIZE)
        #+ geom_text_repel(color=COLOR_SAMPLES_DEF, aes(label= rownames(df)), size = PCH_TEXT_SIZE, force=2)
    }else{
      p = NULL
    }
  }
  p = plotSpace(rgcca, df, "Samples", response, "Response", compX, compY, i_block, p)

  # remove legend if missing
  if (is.null(group)){
    p + theme(legend.position = "none")
  }else
    p
}

getBlocsVariables = function(){
  # Get a vector of block name for each corresponding variable

  rep( names(blocks)[-length(blocks)],
       sapply(blocks[1:(length(blocks)-1)], NCOL))
}

plotVariablesSpace = function(rgcca, compX, compY, i_block=NULL){
  # Projectes a correlation circle in a bi-dimensional space between the coordinates of each variables and its initial value
  # compX: component used for the x-axis
  # compY: component used for the y-axis
  # i_block : index of the block

  if ( is.null(i_block) )
    i_block = length(blocks)

  df =  data.frame(
    #correlation matrix with superblock for each variables and each component selected
    sapply ( c(compX:compY), function(x) cor( blocks[[i_block]], rgcca$Y[[i_block]][, x] ) ) ,
    row.names = colnames(blocks[[i_block]])
  )

  # if superblock is selected, color by blocks
  if (  SUPERBLOCK & ( i_block == length(blocks)) )
    color = getBlocsVariables()
  else
    color = rep(1, NROW(df))

  df = data.frame(df, color )

  p = plotSpace(rgcca, df, "Variables", color, "Blocks", compX, compY, i_block) +
    geom_path(aes(x, y), data = circleFun(), col = "grey", size = 1) +
    geom_path(aes(x, y), data = circleFun()/2, col = "grey", size = 1, lty = 2)

  # remove legend if not on superblock
  if (  !SUPERBLOCK || !( i_block == length(blocks) ) )
      p + theme(legend.position = "none")
    else
      p

}

plotSpace = function (rgcca, df, title, group, name_group, compX, compY, i_block, p=NULL){
  # Projectes coordinates of points in a bi-dimensional space
  # df : dataframe
  # title : type of space (variables or samples)
  # group : color the points with a vector
  # name_group : type of groups (Blocs/Response)
  # compX: component used for the x-axis
  # compY: component used for the y-axis
  # i_block : index of the block

  #if (compX > NB_COMP) compX = 1
  #if (compY > NB_COMP) compY = 2

  p = ggplot(df, aes(df[,compX], df[,compY], colour = group)) +
    geom_text(aes(label = rownames(df)), size = PCH_TEXT_SIZE)
    #geom_text_repel(aes(label= rownames(df)), size = PCH_TEXT_SIZE, force=2)

  p + theme_classic() +
    geom_vline(xintercept = 0, col="grey", linetype="dashed", size=1) +
    geom_hline(yintercept = 0, col="grey", linetype="dashed", size=1) +
    labs ( title = paste(title, "space"),
           x = printAxis(rgcca, compX, i_block),
           y = printAxis(rgcca, compY, i_block),
           color = name_group) +
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

plot_biomarkers = function(rgcca, i_comp, n_mark, i_block=NULL){
  # Histogram plot of the n best biomarkers (according to rgcca$a) on the ieme block
  # i_comp : index of the component
  # n_mark : number of best biomarkers to select
  # i_block : index of the block

  # if no specific block is selected, by default, the superblock is selected (or the last one)
  if ( is.null(i_block) )
    i_block = length(blocks)

  # select the weights
  df = rgcca$a[[i_block]]
  # order by decreasing

  if (  SUPERBLOCK & ( i_block == length(blocks) ) )
    df = data.frame(df, color = getBlocsVariables() )

  df = data.frame(df[order(abs(df[,i_comp]), decreasing = TRUE),], order = nrow(df):1)

  # if superblock is selected, color the bar according to their belonging to each blocks
  #TODO: change this with a booleean with/without superblock
  if (  SUPERBLOCK & ( i_block == length(blocks) ) ){
    # color for the text axis
    color2 = df$color; levels(color2) = hue_pal()(length(blocks)-1)
  }else{
    color2 = "black"
  }

  # max threshold for n
  if(NROW(df) >= n_mark) df = df[1:n_mark,]

  if (  SUPERBLOCK & i_block == length(blocks) ){
    p = ggplot(df, aes(order, df[,i_comp], fill = color))
  }else{
    p = ggplot(df, aes(order, df[,i_comp]))
  }

    p = plotHistogram(p, df, "Variable weights", as.character(color2))
    p + labs (subtitle = printAxis(rgcca, i_comp, i_block))
}

plotAVE = function(rgcca, i_comp){
  # Histogram plot of the Average Variance Explained for each blocks ordered deacreasingly
  # i_comp : index of the component

  df = Reduce(rbind, rgcca$AVE$AVE_X)
  rownames(df) = names(blocks)

  # order by decreasing
  df = data.frame(df[order(abs(df[,i_comp]), decreasing = TRUE),], order = nrow(df):1)

  p = ggplot(df, aes(order, df[,i_comp]))
  plotHistogram(p, df, "Average Variance Explained", "black")
}

plotHistogram = function(p, df, title, color){
  # Default font for a vertical barplot
  # p: ggplot object
  # df : dataframe
  # title: graphic title
  # color: vector of character corresponding to colors for the rows

    p +
    #TODO: if NB_ROW > X, uncomment this
    #geom_hline(yintercept = c(-.5,.5), col="grey", linetype="dotted", size=1) +
    geom_hline(yintercept = 0, col="grey", size=1) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_x_continuous(breaks = df$order, labels = rownames(df)) +
    labs(
      title = title,
      x = "", y = "",
      fill = "Blocks") +
    theme_classic() +
    theme_perso() +
    theme(
      axis.text.y = element_text(size = AXIS_TEXT_SIZE, face = AXIS_FONT, color = color),
      axis.text.x = element_text(size = AXIS_TEXT_SIZE, face = AXIS_FONT, color = "darkgrey"),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, size = 16, face = "italic"))
}
