#' Get the components of the analysis
#' 
#' @inheritParams plotSamplesSpace
#' @inheritParams getVar
#' @param i_block_z An integer giving the index of a list of blocks (another
#' one, different from the one used in i_block)
#' @return A matrix containg each selected components and an associated response
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca.res = rgcca.analyze(blocks)
#' response = factor( apply(Russett[, 9:11], 1, which.max),
#'                   labels = colnames(Russett)[9:11] )
#' getComponents(rgcca.res, as.matrix(response))
#' response = as.matrix(runif(nrow(blocks[[1]])))
#' row.names(response) = row.names(blocks[[1]])
#' getComponents(rgcca.res, response)
#' @export
getComponents <- function(
    rgcca,
    resp = rep(1, nrow(rgcca$Y[[1]])),
    comp_x = 1,
    comp_y = 2,
    comp_z = NULL,
    i_block = length(rgcca$Y),
    i_block_y = i_block,
    i_block_z = i_block,
    predicted = NULL){
    
    df <- data.frame(
        rgcca$Y[[i_block]][, comp_x],
        rgcca$Y[[i_block_y]][, comp_y],
        rgcca$Y[[i_block_z]][, comp_z]
    )
    
    resp <- as.matrix(resp)

    if (!is.null(predicted)) {

            df <- rbind(df, predicted[[2]][[i_block]][, c(comp_x, comp_y, comp_z)])
            resp <- rep(c("obs", "pred"), each = nrow(rgcca$Y[[1]]))

    } else if (length(unique(as.matrix(resp))) > 1) {
        names <- row.names(resp)
        resp <- apply(as.matrix(resp), 1, as.character)

        if (!is.null(names)) {

            resp <- as.matrix(resp, row.names = names)
            name_blocks <- row.names(rgcca$Y[[i_block]])
            diff_column <- setdiff(name_blocks, names)

            if (identical(diff_column, name_blocks)) {
                warning("No match has been found with the row names of the group file.")
                resp <- rep("NA", nrow(df))

            } else {
                if (length(diff_column) > 0) {
                    resp[diff_column] <- "NA"
                    names(resp)[names(resp) == ""] <- names
                } else {
                    names(resp) <- names
                }
                resp <- resp[row.names(rgcca$Y[[i_block]])]
            }
        } else {
            warning("No row names have been found in the group file.")
            resp <- rep("NA", nrow(df))
        }
    }

    if ((!unique(isCharacter(as.vector(resp))) &&
        length(unique(resp)) > 5) || 
            unique(resp) == 1 ) {
        resp[resp == "NA"] <- NA
        df$resp <- as.numeric(resp)
    }else
        df$resp <- as.factor(resp)

    return(df)
}
