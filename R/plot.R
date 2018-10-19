
circleFun = function(center = c(0, 0), diameter = 2, npoints = 100) {
  # creates x,y coordinates for a circle
  
  r = diameter/2
  tt = seq(0, 2 * pi, length.out = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[2] + r * sin(tt)
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
  
  if(unique(isCharacter(as.vector(response)))){
    p = ggplot(df, aes(df[,1], df[,2], colour = response)) +
      geom_text_repel(aes(label= rownames(df)), size = 3, force=2)
  } 
  else{
    p = ggplot(df, aes(df[,1], df[,2], alpha=(response - min(response)) / max(response - min(response)))) +
      scale_alpha_continuous(
        name = "Response",
        breaks = seq(0,1,.25),
        labels = round(quantile(response))
      ) + geom_text_repel(color=COLOR_SAMPLES_DEF, aes(label= rownames(df)), size = 3, force=2)
  }
  
  p + theme_classic() +
    geom_vline(xintercept = 0, col="grey", linetype="dashed", size=1) + 
    geom_hline(yintercept = 0, col="grey", linetype="dashed", size=1) + 
    labs ( title = paste(title, "space"),
           x = printAxis(comp1), 
           y = printAxis(comp2),
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

plot_biomarkers = function(df, comp, n){
  
  df = data.frame(df[order(abs(df[,comp]), decreasing = TRUE),], order = nrow(df):1)
  color2=df$color; levels(color2)=hue_pal()(length(blocks)-1)
  if(NROW(df) >= n) df = df[1:n,]
  
  ggplot(df, aes(order, df[,comp], fill = color)) +
    #geom_hline(yintercept = c(-.5,.5), col="grey", linetype="dotted", size=1) + 
    geom_hline(yintercept = 0, col="grey", size=1) +
    geom_bar(stat = "identity") +
    coord_flip() + 
    scale_x_continuous(breaks=df$order, labels=rownames(df)) +
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
