#' Compute bootstrap (internal)
#'
#' Internal function for computing boostrap of RGCCA
#'
#' @inheritParams rgcca.analyze
#' @inheritParams plot_var_2D
#' @return A list of RGCCA bootstrap weights
#' @examples
#' library(RGCCA)
#' data("Russett")
#' blocks = list(agriculture = Russett[, seq(3)], industry = Russett[, 4:5],
#'     politic = Russett[, 6:11] )
#' rgcca_out = rgcca.analyze(blocks)
#' bootstrap_k(rgcca_out)
#' bootstrap_k(rgcca_out, lapply(blocks, scale), superblock = FALSE)
#' @export
bootstrap_k <- function(
    rgcca,
    A = NULL,
    C = 1 - diag(length(A)),
    tau = rep(1, length(A)),
    ncomp = rep(2, length(A)),
    scheme = "factorial",
    init = "svd",
    bias = TRUE,
    tol = 1e-08,
    type = "rgcca",
    superblock = TRUE) {

    if (is.null(A))
        blocks.all <- rgcca$blocks
    else
        blocks.all <- A

    if (is.null(A)) {
        A <- rgcca$blocks
        C <- rgcca$blocks
        ncomp <- rgcca$ncomp
        scheme <- rgcca$scheme
        bias <- rgcca$bias
        tol <- rgcca$tol
        superblock <- rgcca$superblock
        type <- class(rgcca)

        if (is(rgcca, "sgcca"))
            tau <- rgcca$c1
        else
            tau <- rgcca$tau

        if (superblock) {
            A <- A[-length(A)]
            C <- NULL
        }
    }

    # Shuffle rows
    id_boot <- sample(NROW(A[[1]]), replace = TRUE)

    if (any(sapply(A, function(x) is.null(attr(x, 'scaled:center')))))
            stop("Blocks should be scaled before performing bootstraps.")
    else
        boot_blocks <- lapply(A, function(x) scale2(x[id_boot, ], scale = FALSE))

    boot_blocks <- remove_null_sd(boot_blocks)

    # Get boostraped weights
    w <- rgcca.analyze(
        boot_blocks,
        C,
        superblock = superblock,
        tau = tau,
        ncomp = ncomp,
        scheme = scheme,
        scale = FALSE,
        type = type,
        verbose = FALSE,
        init = init,
        bias = bias,
        tol = tol
    )$a

    # Add removed variables
    missing_var <- lapply(
        seq(length(w)),
        function(x) setdiff(colnames(blocks.all[[x]]), rownames(w[[x]]))
    )

    missing_tab <- lapply(
        seq(length(w)),
        function(x)
            matrix(
                0,
                length(missing_var[[x]]),
                ncomp[x],
                dimnames = list(missing_var[[x]], seq(ncomp[x]))
        ))

    # bug mapply with pca
    w <- lapply(seq(length(missing_tab)), function(x) {
        if (NROW(missing_tab[[x]]) != 0)
            rbind(w[[x]], missing_tab[[x]])
        else
            w[[x]]
        })

    w <- lapply(seq(length(w)), function(x) w[[x]][colnames(blocks.all[[x]]), ])

    names(w) <- names(blocks.all)

    w <- check_sign_comp(rgcca, w)

    return(w)
}
