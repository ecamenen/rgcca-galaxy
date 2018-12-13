
getFileName = function(fi) {
  # get prefix part from a file
  fo = unlist(strsplit(fi, "/"))
  fo = fo[length(fo)]
  fo = unlist(strsplit(fo, "[.]"))[1]
}

loadData = function(fi, fo = fi, row.names = NULL, h = F) {
  # create a dataset object from a file loading fi: input file name fo:
  # dataset object name
  data = as.matrix(read.table(fi, sep = SEPARATOR, h = h, row.names = row.names, na.strings = "NA"))
  assign(fo, data, .GlobalEnv)
  # TODO: catch warning missing \n at the end of the file
}

loadExcel = function(fi, fo = fi, row.names = NULL, h = F) {
  data = read.xlsx2(opt$datasets, fi, header = h)
  checkQuantitative(data[, -row.names], opt$datasets)
  data2 = as.matrix(as.data.frame(lapply(data[-row.names], function(x) as.numeric(as.vector(x)))))
  row.names(data2) = data[, row.names]
  assign(fo, data2, .GlobalEnv)
}

save = function(f, p) {
  pdf(f, width = 10, height = 8)
  plot(p)
  suprLog = dev.off()
}

parseList = function(l) {
  # l: string of characters separated by , out: vector remove white space
  l = gsub(" ", "", l)
  # split by ,
  unlist(strsplit(l, ","))
}

checkQuantitative = function(df, fo) {
  QUALITATIVE = unique(unique(isCharacter(as.matrix(df))))
  if (length(QUALITATIVE) > 1 || QUALITATIVE) {
    msg = paste(fo, "file contains qualitative data. Please, transform them in a disjunctive table.")
    if (!HEADER)
      msg = paste(msg, MSG_HEADER, sep = "")
    stop(paste(msg, "\n"), call. = FALSE)
  }
}

setBlocks = function() {
  # Create a list object of blocks from files loading
  # Output: list of dataframe (blocks)

  # Parse args containing files path
  isXls <- (length(grep("xlsx?", opt$datasets)) == 1)
  # test if extension filename is xls
  if (!isXls) {
    # if it is not, parse the name of file from the arg list
    blocksFilename = parseList(opt$datasets)
  } else {
    # if xls, check file exists
    checkFile(opt$datasets)
    # load the xls
    wb = loadWorkbook(opt$datasets)
    # load the blocks
    blocksFilename = names(getSheets(wb))
  }

  # Parse optional names of blocks
  if (!is.null(opt$names))
    # default name is filename, otherwise, the user could name the blocs
    blocksName = parseList(opt$names)

  # Load each dataset
  blocks = list()
  for (i in 1:length(blocksFilename)) {

    if (!isXls) {
      # if not an xls, file exist test is done here
      fi = blocksFilename[i]
      checkFile(fi)
    }

    #Get names of blocs
    if (!is.null(opt$names))
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
      loadData(fi, fo, 1, HEADER)
    else
      loadExcel(blocksFilename[i], fo, 1, HEADER)

    #if one-column file, it is a tabulation error
    if (NCOL(get(fo)) == 0)
      stop(paste(fo, "block file has an only-column. Check the --separator [by default: 1 for tabulation].\n"),
           call. = FALSE)

    checkQuantitative(get(fo), fo)

    blocks[[fo]] = get(fo)
  }

  if (length(unique(sapply(1:length(blocks), function(x) NROW(blocks[[x]])))) > 1)
    stop("The number of rows is different among the blocks.\n", call. = FALSE)
  #print(names(blocks[[3]])[99])
  #blocks[[3]] = blocks[[3]][, -99]

  if( SUPERBLOCK )
    blocks[["Superblock"]] = Reduce(cbind, blocks)
  #blocks[["Superblock"]] = blocks[["Superblock"]][, -242]

  return(blocks)
}

checkConnection = function(c) {
  # cm: connection matrix unname to avoid taking account the automatic
  # names of column
  if (!isSymmetric.matrix(unname(c)))
    stop("The connection file must be a symmetric matrix.\n", call. = FALSE)
  n = length(blocks)
  if (NCOL(c) != n)
    stop(paste("The number of rows/columns of the connection matrix file must be equals to the number of files in the dataset + 1 (",
               n, ").\n", sep = ""), call. = FALSE)
  d = unique(diag(c))
  if (length(d) != 1 || d != 0)
    stop("The diagonal of the connection matrix file must be 0.\n", call. = FALSE)
  x = unique(c %in% c(0, 1))
  if (length(x) != 1 || x != T)
    stop("The connection file must contains only 0 or 1.\n", call. = FALSE)

}

setConnection = function() {
  # default settings of connection_matrix matrix
  if (is.null(opt$connection)) {
    seq = 1:(length(blocks) - 1)
    connection_matrix = matrix(0, length(blocks), length(blocks))
    connection_matrix[length(blocks), seq] <- connection_matrix[seq, length(blocks)] <- 1
  } else {
    loadData(opt$connection, "connection_matrix", h = F)
  }
  checkConnection(connection_matrix)
  return(connection_matrix)
}

setResponse = function() {
  # create a dataset object from a file loading containg the response
  if ("response" %in% names(opt)) {
    loadData(opt$response, "response", 1, HEADER)
    if (NROW(blocks[[1]]) != NROW(response)) {
      msg = "The number of rows of the response file is different from those of the blocks."
      if (HEADER)
        msg = paste(msg, MSG_HEADER, sep = "")
      stop(paste(msg, "\n"), call. = FALSE)
    }
    QUALITATIVE = unique(isCharacter(response))
    if (length(QUALITATIVE) > 1)
      stop("Please, select a response file with either qualitative data only or quantitative data only. The header must be disabled for quantitative data and activated for disjunctive table.\n",
           call. = FALSE)
    if (NCOL(response) > 1) {
      DISJUNCTIVE = unique(apply(response, 1, sum))
      if (length(DISJUNCTIVE) == 1 && unique(response %in% c(0, 1)) && DISJUNCTIVE == 1) {
        response2 = factor(apply(response, 1, which.max))
        if (HEADER) {
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

isCharacter = function(df) {
  options(warn = -1)
  # is. character() consider a string with '1.2' as a character, not this
  # function NA are produced by converting a character into an integer
  # as.vector, avoid factors of character in integer without NA

  # NA tolerance :
  # df = na.omit(df)
  if (is.matrix(df))
    test = sapply(1:NCOL(df), function(x) unique(is.na(as.integer(as.vector(df[, x])))))
  else
    test = unique(is.na(as.integer(as.vector(df))))

  options(warn = 0)
  return(test)
}
