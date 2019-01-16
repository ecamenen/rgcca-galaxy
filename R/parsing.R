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
MSG_HEADER = " Possible mistake: header parameter is disabled, check if the file doesn't have one."
ROW_NAMES = 1 # column of row names

#' File name from a path
#'
#' Get the file name from a path
#'
#' @param fi A character giving the path of a file
#' @return A character
#' @examples
#' fi = "/name.lastname/dirPath/fileName.tsv"
#' getFileName(fi)
#' # fileName
#' @export getFileName
getFileName = function(fi) {

  fo = unlist(strsplit(fi, "/"))
  fo = fo[length(fo)]
  unlist(strsplit(fo, "[.]"))[1]
}

#' Creates a matrix from loading a file
#'
#' @param f A character giving the file name
#' @param sep A character giving the column separator
#' @param row.names A vector of characters giving the names of the rows
#' @param h A bolean giving the presence or the absence of the header
#' @return A matrix containing the loaded file
#' @examples
#' \dontrun{
#' loadData("data/agriculture.tsv")
#' }
#' @export loadData
loadData = function(f, sep = "\t", row.names = 1, h = TRUE) {

  as.matrix(read.table(f, sep = sep, header = h, row.names = row.names, na.strings = "NA"))
  # TODO: catch warning missing \n at the end of the file
}

#' Creates a data frame from an Excel file loading
#'
#' @param f A character giving the file name
#' @param sheet A character giving the sheet name
#' @param row.names A vector of characters giving the names of the rows
#' @param h A bolean giving the presence or the absence of the header
#' @return A matrix containing the loaded file
#' @examples
#' \dontrun{
#' loadExcel("data/blocks.xlsx", "industry")
#' }
#' @export loadExcel
loadExcel = function(f, sheet, row.names = 1, h = TRUE) {
  # TODO: opt$dataset in arg

  df = read.xlsx2(f, sheet, header = h)
  checkQuantitative(df[, -row.names], opt$datasets, h)
  df2 = as.matrix(as.data.frame(lapply(df[-row.names], function(x) as.numeric(as.vector(x)))))
  row.names(df2) = df[, row.names]

  return (df2)
}

#' Save a ggplot object
#'
#' Save a ggplot in various output formats
#'
#' @param f A character giving the name of a file
#' @param p A ggplot object
#' @examples
#' library("ggplot2")
#' df = as.data.frame(matrix(runif(20), 10, 2))
#' p = ggplot(df, aes(df[, 1], df[, 2]))
#' savePlot("Rplot.png", p)
#' @export savePlot
savePlot = function(f, p) {

  # get suffixe of filename
  format = unlist(strsplit(f, '.', fixed="T"))
  format = format[length(format)]

  # dynamic loading of function depending of the extension
  if (format == "dat")
    func = pdf
  else
    func = get(format)

  # save
  if (format %in% c("pdf", "dat") ) func(f, width = 10, height = 8)
  else func(f, width = 10, height = 8, units = "in", res = 200)

  plot(p)
  suprLog = dev.off()
}

#' Convert a character in a vector
#'
#' @param s A character separated by comma
#' @return A vector of characters whitout spaces
#' @examples
#' s = "1,2, 3"
#' parseList(s)
#' @export parseList
parseList = function(s) {

  s = gsub(" ", "", s)
  # split by comma
  unlist(strsplit(s, ","))
}

#' Check if a dataframe contains no qualitative variables
#'
#' @param df A dataframe or a matrix
#' @param fo A character giving the name of the tested file
#' @param h A bolean giving either the presence (TRUE) or absence (FALSE) of a header
#' @examples
#' df = matrix(runif(20), 10, 2)
#' checkQuantitative(df, "data")
#' \dontrun{
#' df[,2] = LETTERS[1:10]
#' checkQuantitative(df, "data", TRUE)
#' # Error
#' }
#' @export checkQuantitative
checkQuantitative = function(df, fo, h = FALSE) {
  qualitative = unique(unique(isCharacter(as.matrix(df))))
  if (length(qualitative) > 1 || qualitative) {
    msg = paste(fo, "file contains qualitative data. Please, transform them in a disjunctive table.")
    if (!h)
      msg = paste(msg, MSG_HEADER, sep = "")
    stop(paste(msg, "\n"), call. = FALSE)
  }
}

checkFile = function (f){
  # Check the existence of a path
  # f: A character giving the path of a file

  if(!file.exists(f)){
    stop(paste(f, " file does not exist\n", sep=""), call.=FALSE)
  }
}

#' Create a list of matrix from loading files corresponding to blocks
#'
#' @param superblock A boolean giving the presence (TRUE) / absence (FALSE) of a superblock
#' @param file A character giving the path of a file used as a response
#' @param names A character giving a list of names for the blocks
#' @param sep A character giving the column separator
#' @param header A bolean giving the presence or the absence of the header
#' @return A list matrix corresponding to the blocks
#' @examples
#' \dontrun{
#' setBlocks (TRUE, "data/agriculture.tsv,data/industry.tsv,data/politic.tsv", "agric,ind,polit")
#' }
#' @export setBlocks
setBlocks = function(superblock, file, names = NULL, sep = "\t", header = TRUE) {

  # Parse args containing files path
  isXls <- (length(grep("xlsx?", file)) == 1)
  # test if extension filename is xls
  if (!isXls)
    # if it is not, parse the name of file from the arg list
    blocksFilename = parseList(file)
  else {
    # if xls, check file exists
    checkFile(file)
    # load the xls
    wb = loadWorkbook(file)
    # load the blocks
    blocksFilename = names(getSheets(wb))
  }

  # Parse optional names of blocks
  if (!is.null(names))
    # default name is filename, otherwise, the user could name the blocs
    blocksName = parseList(names)

  # Load each dataset
  blocks = list()
  for (i in 1:length(blocksFilename)) {

    if (!isXls) {
      # if not an xls, file exist test is done here
      fi = blocksFilename[i]
      checkFile(fi)
    }

    #Get names of blocs
    if (!is.null(names))
      # names of blocks are those parsed from args
      fo = getFileName(blocksName[i])
    else {
      if (!isXls)
        # if not xls, the name is the files without the extension .tsv
        fo = getFileName(fi)
      else
        # for xls, the names are those of the sheets
        fo = blocksFilename[i]
    }

    #load the data
    if (!isXls)
      df = loadData(fi, sep, ROW_NAMES, header)
    else
      df = loadExcel(file, blocksFilename[i], ROW_NAMES, header)

    #if one-column file, it is a tabulation error
    if (NCOL(df) == 0)
      stop(paste(fo, "block file has an only-column. Check the separator [by default: tabulation].\n"),
           call. = FALSE)

    checkQuantitative(df, fo, header)

    blocks[[fo]] = df
  }

  if (length(unique(sapply(1:length(blocks), function(x) NROW(blocks[[x]])))) > 1)
    stop("The number of rows is different among the blocks.\n", call. = FALSE)
  #print(names(blocks[[3]])[99])
  #blocks[[3]] = blocks[[3]][, -99]

  if( superblock )
    blocks[["Superblock"]] = Reduce(cbind, blocks)
  #blocks[["Superblock"]] = blocks[["Superblock"]][, -242]

  return(blocks)
}

#' Check the format of the connection matrix
#'
#' @param c A symmetric matrix containing 1 and 0
#' @param blocks A list of matrix
#' @export checkConnection
checkConnection = function(c, blocks) {

  if (!isSymmetric.matrix(unname(c)))
    stop("The connection file must be a symmetric matrix.\n", call. = FALSE)
  n = length(blocks)
  if (NCOL(c) != n)
    stop(paste("The number of rows/columns of the connection matrix file must be equals to ",
               n,
               " (the number of blocks in the dataset, +1 with a superblock by default).\n", sep = ""),
         call. = FALSE)
  d = unique(diag(c))
  if (length(d) != 1 || d != 0)
    stop("The diagonal of the connection matrix file must be 0.\n", call. = FALSE)
  x = unique(c %in% c(0, 1))
  if (length(x) != 1 || x != T)
    stop("The connection file must contains only 0 or 1.\n", call. = FALSE)

}

#' Create a matrix from loading a file corresponding to a connection between the blocks
#'
#' @param blocks A list of matrix
#' @param file A character giving the path of a file used as a response
#' @param sep A character giving the column separator
#' @return A matrix corresponding to the response
#' @examples
#' \dontrun{
#' blocks = lapply(1:4, function(x) matrix(runif(47 * 5), 47, 5))
#' setConnection (blocks, "data/connection.tsv")
#' }
#' @export setConnection
setConnection = function(blocks, file = NULL, sep = "\t") {

  if (is.null(file)) {

    seq = 1:(length(blocks) - 1)
    connection = matrix(0, length(blocks), length(blocks))
    connection[length(blocks), seq] <- connection[seq, length(blocks)] <- 1

  } else {
    connection = loadData(file, sep, NULL, FALSE)
  }

  checkConnection(connection, blocks)

  return(connection)
}


#' Create a matrix from loading a file corresponding to the response
#'
#' @param blocks A list of matrix
#' @param file A character giving the path of a file used as a response
#' @param sep A character giving the column separator
#' @param header A bolean giving the presence or the absence of the header
#' @return A matrix corresponding to the response
#' @examples
#' \dontrun{
#' blocks = lapply(1:3, function(x) matrix(runif(47 * 5), 47, 5))
#' setResponse (blocks, "data/response3.tsv")
#' }
#' @export setResponse
setResponse = function(blocks, file = NULL, sep = "\t", header = TRUE) {

  if (!is.null(file)) {
    response = loadData(file, sep, ROW_NAMES, header)

    if (NROW(blocks[[1]]) != NROW(response)) {
      msg = paste("The number of rows of the response file (", NROW(response), ") is different from those of the blocks (", NROW(blocks[[1]]), ").", sep="")
      if (header)
        msg = paste(msg, MSG_HEADER, sep = "")
      stop(paste(msg, "\n"), call. = FALSE)
    }

    qualitative = unique(isCharacter(response))

    if (length(qualitative) > 1)
      stop("Please, select a response file with either qualitative data only or quantitative data only. The header must be disabled for quantitative data and activated for disjunctive table.\n",
           call. = FALSE)

    if (NCOL(response) > 1) {
      disjunctive = unique(apply(response, 1, sum))

      if (length(disjunctive) == 1 && unique(response %in% c(0, 1)) && disjunctive == 1) {
        response2 = factor(apply(response, 1, which.max))
        if (header) {
          levels(response2) = colnames(response)
        }
        response = as.character(response2)

      } else {
        response = response[, 1]
        warning("There is multiple columns in the response file. By default, only the first one is taken in account.\n",
                call. = FALSE)
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
#' @return A bolean for the presence (FALSE) or the absence (TRUE) of at least one quantitative variable
#' @examples
#' x = matrix(c(runif(10), LETTERS[1:10]), 10, 2)
#' isCharacter(x)
#' # FALSE TRUE
#' isCharacter(LETTERS[1:10])
#' # TRUE
#' @export isCharacter
isCharacter = function(x) {

  options(warn = -1)
  # is. character() consider a string with '1.2' as a character, not this
  # function NA are produced by converting a character into an integer
  # as.vector, avoid factors of character in integer without NA

  # NA tolerance :
  # x = na.omit(x)
  if (is.matrix(x))
    test = sapply(1:NCOL(x), function(i) unique(is.na(as.integer(as.vector(x[, i])))))
  else
    test = unique(is.na(as.integer(as.vector(x))))

  options(warn = 0)
  return(test)
}
