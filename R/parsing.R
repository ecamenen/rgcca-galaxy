# Author: Etienne CAMENEN
# Date: 2019
# Contact: arthur.tenenhaus@l2s.centralesupelec.fr
# Key-words: omics, RGCCA, multi-block
# EDAM operation: analysis, correlation, visualisation
#
# Abstract: Performs multi-variate analysis (PCA, CCA, PLS, R/SGCCA, etc.)
# and produces textual and graphical outputs (e.g. variables and individuals
# plots).

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

warning <- function(message, call = sys.call(-1)) {
    base::warning(message, call. = FALSE, immediate. = TRUE)
}

stop <- function(message,
                exit_code = "1",
                call = NULL) {
    base::stop(structure(
        class = c(exit_code, "simpleError", "error", "condition"),
        list(message = message, call. = NULL)
    ))
}

# Global settings
MSG_HEADER <-
" Possible mistake: header parameter is disabled, check if the file doesn't have one."
ROW_NAMES <- 1  # column of row names

#' File name from a path
#'
#' Get the file name from a path
#'
#' @param fi A character giving the path of a file
#' @return A character for the name of the file
#' @examples
#' fi = '/name.lastname/dirPath/fileName.tsv'
#' getFileName(fi)
#' # fileName
#' @export
getFileName <- function(fi) {
    if (!is.null(fi)) {
        fo <- unlist(strsplit(fi, "/"))
        fo <- fo[length(fo)]
        unlist(strsplit(fo, "[.]"))[1]
    }
}

# Print warning if file size over
checkFileSize <- function(filename) {
    size <- file.size(filename)
    if (size > 5e+06)
        # warning(paste0('The size of ', filename, ' is over 5 Mo (',
        #  round(size / 1E6, 1), ' Mo). File loading could take some times...'),
        message("File loading in progress ...")
}

convertMatrixNumeric <- function(df) {
    matrix(
        sapply(seq_len(nrow(df) * ncol(df)), function(i)
            tryCatch({
                as.numeric(df[i])
            }, warning = function(e)
                NA)),
        nrow(df),
        ncol(df),
        dimnames = list(row.names(df),
                        colnames(df))
    )
}

#' Creates a matrix from loading a file
#'
#' @param f A character giving the file name
#' @param sep A character giving the column separator
#' @param rownames An integer corresponding to the column number of the
#' row names (NULL otherwise)
#' @param h A bolean giving the presence or the absence of the header
#' @return A matrix containing the loaded file
#' @examples
#' \dontrun{
#' loadData('data/agriculture.tsv')
#' }
#' @export
loadData <- function(f, sep = "\t", rownames = 1, h = TRUE) {

    if (!is.null(rownames) && rownames < 1)
        rownames <- NULL

    func <- function(x = rownames)
            as.matrix(read.table(
                f,
                sep = sep,
                header = h,
                row.names = x,
                na.strings = "NA",
                dec = ","
            ))

    tryCatch({
        func()
    }, error = function(e) {
        if (e$message == "duplicate 'row.names' are not allowed")
            func(NULL)
    })

}

# Creates a data frame from an Excel file loading
#
# @param f A character giving the file name
# @param sheet A character giving the sheet name
# @param rownames An integer corresponding to the column number of the row
#' names (NULL otherwise)
# @param h A bolean giving the presence or the absence of the header
# @param num A bolean giving the presence or the absence of numerical values
# @return A matrix containing the loaded file
# @examples
# \dontrun{
# loadExcel("data/blocks.xlsx", "industry")
# }
# @export loadExcel
# loadExcel = function(f, sheet, rownames = 1, h = TRUE, num = TRUE) {
#
#   if (!is.null(rownames) && rownames < 1)
#     rownames = NULL
#
#   df = read.xlsx(f, sheet, header = h, startRow = 1)
#
#   if(!is.null(rownames)){
#     names = df[, rownames]
#     df = df[, -rownames]
#   }
#
#   if(num)
#     df = as.data.frame(lapply(df, function(x) as.numeric(as.vector(x))))
#
#   df = as.matrix(df)
#
#   if(!is.null(rownames))
#     row.names(df) = names
#
#   return (df)
# }

#' Save a ggplot object
#'
#' Save a ggplot in various output formats
#'
#' @param f A character giving the name of a file
#' @param p A ggplot object
#' @examples
#' library('ggplot2')
#' df = as.data.frame(matrix(runif(20), 10, 2))
#' p = ggplot(df, aes(df[, 1], df[, 2]))
#' #savePlot('Rplot.png', p)
#' @export
savePlot <- function(f, p) {

    # get suffixe of filename
    format <- unlist(strsplit(f, ".", fixed = "TRUE"))
    format <- format[length(format)]

    # dynamic loading of formattion depending of the extension
    if (format == "dat")
        formatFunc <- pdf
    else
        formatFunc <- get(format)

    # save
    if (format %in% c("pdf", "dat"))
        formatFunc(f, width = 10, height = 8)
    else
        formatFunc(
            f,
            width = 10,
            height = 8,
            units = "in",
            res = 200
        )

    if (is.function(p))
        p()
    else
        plot(p)

    invisible(dev.off())
}

#' Convert a character in a vector
#'
#' @param s A character separated by comma
#' @return A vector of characters whitout spaces
#' @examples
#' s = '1,2, 3'
#' parseList(s)
#' @export
parseList <- function(s) {
    s <- gsub(" ", "", s)
    # split by comma
    unlist(strsplit(s, ","))
}

#' Check if a dataframe contains no qualitative variables
#'
#' @param df A dataframe or a matrix
#' @param fo A character giving the name of the tested file
#' @param h A bolean giving either the presence (TRUE) or absence (FALSE) of
#'  a header
#' @examples
#' df = matrix(runif(20), 10, 2)
#' checkQuantitative(df, 'data')
#' \dontrun{
#' df[,2] = LETTERS[seq_len(10)]
#' checkQuantitative(df, 'data', TRUE)
#' # Error
#' }
#' @export
checkQuantitative <- function(df, fo, h = FALSE) {
    qualitative <- unique(unique(isCharacter(as.matrix(df))))

    if (length(qualitative) > 1 || qualitative) {
        msg <- paste(
                fo,
                "file contains qualitative data. Please, transform them in a disjunctive table."
            )

        if (!h)
            msg <- paste0(msg, MSG_HEADER)

        stop(paste(msg, "\n"), exit_code = 100)
    }

}

checkFile <- function(f) {
    # Check the existence of a path f: A character giving the path of a file

    if (!file.exists(f))
        stop(paste0(f, " file does not exist."), exit_code = 101)

}

#' Create a list of matrix from loading files corresponding to blocks
#'
#' @param file A character giving the path of a file used as a response
#' @param names A character giving a list of names for the blocks
#' @param sep A character giving the column separator
#' @param header A bolean giving the presence or the absence of the header
#' @param rownames An integer corresponding to the column number of the row
#' names (NULL otherwise)
#' @return A list matrix corresponding to the blocks
#' @examples
#' \dontrun{
#' setBlocks (TRUE,
#'     "data/agriculture.tsv,data/industry.tsv,data/politic.tsv",
#'     "agric,ind,polit")
#' }
#' @export
setBlocks <- function(file,
    names = NULL,
    sep = "\t",
    header = TRUE,
    rownames = ROW_NAMES) {
    
    # Parse args containing files path
    isXls <- (length(grep("xlsx?", file)) == 1)
    # test if extension filename is xls
    if (!isXls)
        # if it is not, parse the name of file from the arg list
        blocksFilename <- parseList(file)
    else {
        # # if xls, check file exists
        # checkFile(file)
        # # load the xls
        # wb = loadWorkbook(file)
        # # load the blocks
        # blocksFilename = names(getSheets(wb))
    }

    # Parse optional names of blocks
    if (!is.null(names))
        # default name is filename, otherwise, the user could name the blocs
        blocksName <- parseList(names)

    # Load each dataset
    blocks <- list()
    for (i in seq_len(length(blocksFilename))) {
        if (!isXls) {
            # if not an xls, file exist test is done here
            fi <- blocksFilename[i]
            checkFile(fi)
        }

        #Get names of blocs
        if (!is.null(names))
            # names of blocks are those parsed from args
            fo <- getFileName(blocksName[i])
        else {
            if (!isXls)
                # if not xls, the name is the files without the extension .tsv
                fo <- getFileName(fi)
            else
                # for xls, the names are those of the sheets
                fo <- blocksFilename[i]
        }

        #load the data
        if (!isXls) {
            checkFileSize(fi)
            df <- loadData(fi, sep, rownames, header)
        }
        # }else{
        #   checkFileSize(file)
        #   df = loadExcel(file, blocksFilename[i], rownames, header)
        # }

        #if one-column file, it is a tabulation error
        if (NCOL(df) == 0)
            stop(paste(fo, "block file has an only-column. Check the separator."),
                exit_code = 102)

        dimnames <- list(row.names(df), colnames(df))
        df <- convertMatrixNumeric(df)

        df <- imputeMean(df)
        
        checkQuantitative(df, fo, header)
        df <- matrix(as.numeric(df), nrow(df), ncol(df), dimnames = dimnames)
        blocks[[fo]] <- df
    }

    nrow <- lapply(blocks, NROW)

    if (length(blocks) > 1)
        blocks <- keepCommonRow(blocks)
    
    blocks <- removeColumnSdNull(blocks)

    for (i in seq_len(length(blocks))) 
        attributes(blocks[[i]])$nrow <- nrow[[i]]

    if (nrow(blocks[[1]]) > 0)
        return(blocks)
    else
        stop("There is no rows in common between the blocks.", exit_code = 108)
}

#' Impute NA
#' 
#' Impute non available by means
#' 
#' @param df A matrix containing non availables data
#' @examples
#' df = cbind(runif(9), runif(9))
#' df = rbind(df, c(NA, NA))
#' imputeMean(df)
#' @return A matrix with imputed values
#' @export
imputeMean <- function(df){
    if (any(is.na(df))) {
        df <- matrix(unlist(
            lapply(seq_len(ncol(df)),
                function(x)
                    unlist(lapply(as.list(df[, x]),
                        function(y)
                            ifelse(is.na(y),
                            mean(unlist(df[, x]), na.rm = TRUE), y))))),
            nrow(df),
            ncol(df), 
            dimnames = list(row.names(df), colnames(df)))
    }
    return(df)
}


#' Check the format of the connection matrix
#'
#' @param c A symmetric matrix containing 1 and 0
#' @param blocks A list of matrix
#' @export
checkConnection <- function(c, blocks) {

    if (!isSymmetric.matrix(unname(c)))
        stop("The connection file must be a symmetric matrix.", exit_code = 103)

    d <- unique(diag(c))
    if (length(d) != 1 || d != 0)
        stop("The diagonal of the connection matrix file must be 0.",
            exit_code = 105)

    x <- unique(c %in% c(0, 1))
    if (length(x) != 1 || x != TRUE)
        stop("The connection file must contains only 0 or 1.", exit_code = 106)

    if (all(c == 0))
        stop("The connection file could not contain only 0.", exit_code = 107)

    n <- length(blocks)
    if (NCOL(c) != n)
        stop(
            paste0(
                "The number of rows/columns of the connection matrix file must be equal to ",
                n,
                " (the number of blocks in the dataset, +1 with a superblock by default)."
            ),
            exit_code = 104
        )

    # TODO: warning if superblock = TRUE

}

#' Create a matrix corresponding to a connection between the blocks
#'
#' @param blocks A list of matrix
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of
#' a superblock
#' @param file A character giving the path of a file used as a response
#' @param sep A character giving the column separator
#' @param rownames An integer corresponding to the column number of the row
#' names (NULL otherwise)
#' @param h A bolean giving the presence or the absence of the header
#' @return A matrix corresponding to the connection between the blocks
#' @examples
#' \dontrun{
#' blocks = lapply(seq_len(4), function(x) matrix(runif(47 * 5), 47, 5))
#' setConnection (blocks, 'data/connection.tsv')
#' }
#' @export
setConnection <- function(blocks,
    superblock = FALSE,
    file = NULL,
    sep = "\t",
    h = FALSE,
    rownames = NULL) {
    

    J <- length(blocks)

    if (superblock) {
        connection <- matrix(0, J, J)
        connection[seq_len(J - 1), J] <- connection[J, seq_len(J - 1)] <- 1

    } else if (is.null(file))
        connection <- 1 - diag(J)
    else {
        isXls <- (length(grep("xlsx?", file)) == 1)

        if (!isXls)
            connection <- loadData(
                    f = file,
                    sep = sep,
                    rownames = rownames,
                    h = h
                )
        # else 
        # connection = loadExcel(f = file, sheet = 1, rownames = rownames, h = h)
    }

    checkConnection(connection, blocks)

    return(connection)
}


#' Create a matrix corresponding to the response
#'
#' @param blocks A list of matrix
#' @param file A character giving the path of a file used as a response
#' @param sep A character giving the column separator
#' @param header A bolean giving the presence or the absence of the header
#' @param rownames An integer corresponding to the column number of the row
#' names (NULL otherwise)
#' @return A matrix corresponding to the response
#' @examples
#' \dontrun{
#' blocks = lapply(seq_len(3), function(x) matrix(runif(47 * 5), 47, 5))
#' setResponse (blocks, 'data/response3.tsv')
#' }
#' @export
setResponse <- function(
    blocks = NULL,
    file = NULL,
    sep = "\t",
    header = TRUE,
    rownames = ROW_NAMES) {
    

    if (!is.null(file)) {
        isXls <- length(grep("xlsx?", file))

        if (!isXls)
            response <- loadData(file, sep, rownames, header)
        # else response = loadExcel(file, 1, rownames, h = header, num = FALSE)

        qualitative <- unique(isCharacter(response))

        if (length(qualitative) > 1)
            stop(
                "Please, select a response file with either qualitative data only or quantitative data only.",
                108
            )

        if (!qualitative)
            response <- convertMatrixNumeric(response)


        if (NCOL(response) > 1) {
            disjunctive <- unique(apply(response, 1, sum))
            
            

            if (length(disjunctive) &&
                unique(disjunctive %in% c(0, 1)) && disjunctive) {
                response2 <- factor(apply(response, 1, which.max))
                if (header) {
                    levels(response2) <- colnames(response)
                }
                response <- as.matrix(data.frame(
                    as.character(response2),
                    row.names = rownames(response)
                    ))

            } else {
                response <- as.matrix(response[, 1])
                warning("There is multiple columns in the response file. By default, only the first one is taken in account.")
            }
        }

        return(response)
    } else {
        return(rep(1, NROW(blocks[[1]])))
    }
}

#' Test for character vector
#'
#' Tests if a dataframe is composed only by qualitative variables
#'
#' @param x A matrix or a vector
#' @return A bolean for the presence (FALSE) or the absence (TRUE) of at least
#' one quantitative variable
#' @examples
#' x = matrix(c(runif(10), LETTERS[seq_len(10)]), 10, 2)
#' isCharacter(x)
#' # FALSE TRUE
#' isCharacter(LETTERS[seq_len(10)])
#' # TRUE
#' @export 
isCharacter <- function(x) {
    # is. character() consider a string with '1.2' as a character, not this function.
    # NA are produced by converting a character into an integer as.vector, avoid
    # factors of character in integer without NA

    # NA tolerance :

    if (is.matrix(x)) {
        test <- sapply(seq_len(NCOL(x)), function(i)
                unique(is.na(
                    tryCatch(
                    as.integer(
                        na.omit(as.vector(x[, i])[as.vector(x[, i]) != "NA"])),
                    warning = function(w)
                        return(NA)
                ))))
    } else
        test <- unique(is.na(
            tryCatch(
                as.integer(na.omit(as.vector(x)[as.vector(x) != "NA"])),
            warning = function(w)
                return(NA)
            )))

    return(test)
}

#' Get the rows with the same names among a list of dataframe
#'
#' @param list_m A list of dataframe
#' @return A vector of character with the common rownames
#' @export
commonRow <- function(list_m) {
    common_row <- row.names(list_m[[1]])
    for (i in 2:length(list_m))
        common_row <- common_row[common_row %in% row.names(list_m[[i]])]
    return(common_row)
}

#' Keep only the rows with the same names among a list of dataframe
#'
#' @param list_m A list of dataframe
#' @return A list of dataframe
#' @export
keepCommonRow <- function(list_m) {

    names <- names(list_m)
    
    common_row <- commonRow(list_m)
    list_m <- lapply(seq_len(length(list_m)), function(x)
            list_m[[x]] <- list_m[[x]][common_row, ])

    names(list_m) <- names
    return(list_m)
}

#' Remove column having a standard deviation equals to 0
#'
#' @param list_m A list of dataframe
#' @return A list of dataframe
#' @examples
#' df = sapply(seq(3), function(x) runif(10))
#' df = cbind(df, rep(1, 10))
#' removeColumnSdNull(list(df))
#' @export
removeColumnSdNull <- function(list_m) {

    names <- names(list_m)

    column_sd_null <- lapply(list_m, function(x)
            which(apply(x, 2, sd) == 0))
    blocks_index <- seq(1, length(list_m))[unlist(lapply(
            column_sd_null, 
            function(x) length(x) > 0))]

    list_m <- lapply(seq_len(length(list_m)), function(x) {
        if (x %in% blocks_index)
            list_m[[x]][, -column_sd_null[[x]]]
        else
            list_m[[x]]
    })

    names(list_m) <- names
    return(list_m)
}

setSuperblock <- function(blocks, superblock = FALSE, type = "rgcca") {

    if (superblock | tolower(type) == "pca") {
        # if(type != 'pca') warnconnection('superblock')
        blocks[["superblock"]] <- Reduce(cbind, blocks)
    }
    return(blocks)
}

setPosPar <- function(opt, blocks, i_resp) {

    J <- length(blocks)
    opt$blocks <- blocks
    opt$block_names <- names(blocks)

    par <- c("blocks", "block_names", "ncomp")
    if (all(opt$tau != "optimal"))
        par[length(par) + 1] <- "tau"

    for (i in seq_len(length(par))) {
        temp <- opt[[par[i]]][[J]]
        opt[[par[i]]][[J]] <- opt[[par[i]]][[i_resp]]
        opt[[par[i]]][[i_resp]] <- temp
    }

    names(opt$blocks) <- opt$block_names

    return(opt)
}


warnConnection <- function(x)
        warning(
            paste0(
                "By using a ",
                x,
                ", all blocks are connected to this block in the connection matrix and the connection file is ignored."
            )
        )

checkSuperblock <- function(opt) {
    if (!is.null(opt$response)) {
        warnConnection("supervized method with a response")
        if (opt$superblock) {
            opt$superblock <- FALSE
            if ("superblock" %in% names(opt))
                warning("In a supervised mode, the superblock corresponds to the response.")
        }
    }
    return(opt)
}
