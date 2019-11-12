#' Plot in 3 dimensions
#' 
#' Plot in 3 dimensions either to visualize the components of an analyse or the variables
#' @inheritParams plotSamplesSpace
#' @inheritParams plotSpace
#' @inheritParams getComponents
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks, ncomp = rep(3, 2))
#' df = getComponents(rgcca.res, comp_z = 3)
#' spacePlot3D(df, i_block = 2)
#' spacePlot3D(df, i_block = 2, text = FALSE)
#' response = factor( apply(Russett[, 9:11], 1, which.max),
#'                   labels = colnames(Russett)[9:11] )
#' response = blocks[[2]][, 1]
#' names(response) = row.names(blocks[[2]])
#' df = getComponents(rgcca.res, response, comp_z = 3)
#' spacePlot3D(df, i_block = 2, text = FALSE)
#' spacePlot3D(df, i_block = 2)
#' df = getVariablesIndexes(rgcca.res, blocks, comp_z = 3, i_block = 1, collapse = TRUE)
#' spacePlot3D(df, i_block = 2, type = "var")
#' @export
spacePlot3D <- function(
    df,
    comp_x = 1,
    comp_y = 2,
    comp_z = 3,
    i_block,
    i_block_y = i_block,
    i_block_z = i_block,
    text = TRUE,
    title = "Sample plot",
    type = "ind") {
    
    if (length(unique(df$resp)) == 1) {
        df$resp = as.factor(rep("a", length(df$resp)))
        midcol = "#cd5b45"
    } else
        midcol = "gray"

    axis <- function(x, i)
        list(
                title = paste0("<i>", printAxis(rgcca.res, x, i), "</i>"),
                titlefont = list(
                        size = AXIS_TITLE_CEX * 0.75
                    )
            )
    
    colorNumeric <- function(x){
        n <- length(x)
         cut(
             x, 
             breaks = n,
             labels = colorRampPalette(c("#A50026", midcol,  "#313695"))(n), 
             include.lowest = TRUE)
    }
    
    subdf <- function(x) 
        df[which(df$resp == levels(df$resp)[x]), ]

    add_trace_manual <- function(p, x){

        l <- levels(df$resp)

        func <- quote(
            add_trace(
                p,
                name = l[x],
                x = ~ subdf(x)[, 1],
                y = ~ subdf(x)[, 2],
                z = ~ subdf(x)[, 3],
                type = "scatter3d",
                showlegend = TRUE
            )
        )
        
        color <- colorGroup(1:length(l))[x]
        
        if (text) {
            func$mode <- "text"
            func$text <- ~row.names(subdf(x))
            func$textfont <- list(
                color = color,
                size = PCH_TEXT_CEX * 4
            )
        }else{
            func$mode <- "markers"
            func$marker <- list(
                color = color,
                size = PCH_TEXT_CEX * 1.5
            )
        }
        
        eval(func)
    }
    
    
    if (!isCharacter(df$resp)) {

        if (text)
            visible <- "legendonly"
        else
            visible <- TRUE

        p <- plot_ly(
            name = "samples",
            x = ~ df[, 1],
            y = ~ df[, 2],
            z = ~ df[, 3],
            mode = "markers",
            type = "scatter3d",
            showlegend = FALSE,
            color = df$resp,
            size = I(200),
            colors = c("#A50026", midcol,  "#313695"),
            visible = visible
        )

        if (text) {
            p <- p %>%
                add_trace(
                    name = "samples",
                    x = ~ df[, 1],
                    y = ~ df[, 2],
                    z = ~ df[, 3],
                    mode = "text",
                    type = "scatter3d",
                    text = ~ row.names(df),
                    textfont = list(
                        color = colorNumeric(df$resp),
                        size = PCH_TEXT_CEX * 4
                    ),
                    showlegend = FALSE,
                    visible = TRUE
                )
        }

    }else{
        p <- plot_ly()
        
        for (i in seq(length(levels(df$resp))))
            p <- p %>% add_trace_manual(i)
    }

    p <- p %>%
        layout(
            autosize = T,
            margin = list(
                l = 50,
                r = 50,
                b = 50,
                t = 100
            ),
            scene = list(
                aspectmode = 'cube',
                xaxis = axis(comp_x, i_block),
                yaxis = axis(comp_y, i_block_y),
                zaxis = axis(comp_z, i_block_z)
            ),
            title = list(
                text = paste0('<b>', title, '</b>'),
                font = list(
                    size = 25 * CEX,
                    face = "bold"
                )
            )
        )

    plot_circle <- function(p, x, y, z){
        df <- cbind(circleFun(), 0)
        add_trace(
            p = p,
            x = df[, x],
            y = df[, y],
            z = df[, z],
            showlegend = FALSE,
            hoverinfo = "none" ,
            mode = "lines",
            type = "scatter3d",
            line = list(color = "grey", width = 4)
        )
    }
        
    if (type == "var")
        p <- p %>% plot_circle(1, 2, 3) %>% plot_circle(1, 3, 2) # %>% plot_circle(3, 2, 1)
    
    return(p)
}
