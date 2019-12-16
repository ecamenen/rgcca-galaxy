# An intern function used by sgcca.permute to perform multiple sgcca with permuted rows
# sgcca.crit(A = blocks, scale = FALSE)
sgcca.crit <- function(
    blocks,
    par = list("ncomps", expand.grid(rep(list(1:2), length(blocks)))),
    connection = 1 - diag(length(blocks)),
    response = NULL,
    tau = rep(1, length(blocks)),
    ncomp = rep(2, length(blocks)),
    scheme = "factorial",
    scale = TRUE,
    init = "svd",
    bias = TRUE,
    tol = 1e-08,
    type = "rgcca",
    superblock = TRUE,
    perm = TRUE,
    n_cores = parallel::detectCores() - 1) {

    if (perm) {
        for (k in 1:length(blocks))
            blocks[[k]] <- blocks[[k]][sample(1:nrow(blocks[[k]])),]
    }

    simplify2array(
        parallel::mclapply(
            seq(NROW(par[[2]])),
            function(i) {
                switch(
                    par[[1]],
                    "ncomp" = {
                        ncomp <- par[[2]][i, ]
                    },
                    "c1" = {
                        tau <- par[[2]][i, ]
                        type <- "sgcca"
                })
                crit <- rgcca.analyze(
                        blocks = blocks,
                        connection = connection,
                        response = response,
                        superblock = superblock,
                        tau = tau,
                        ncomp = ncomp,
                        scheme = scheme,
                        scale = scale,
                        type = type,
                        init = init,
                        bias = bias,
                        tol = tol
                    )$crit
                return(mean(sapply(crit, function(x) x[length(x)])))
            },
        mc.cores = n_cores))
}
