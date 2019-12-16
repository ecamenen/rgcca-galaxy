crossvalidation.gcca <- function(
    rgcca,
    bloc_to_pred = names(rgcca$blocks)[1],
    validation = "loo",
    type = "regression",
    fit = "lm",
    new_scaled = TRUE,
    k = 5,
    rep = 10,
    nb_cores = parallel::detectCores() - 1) {

    match.arg(validation, c("test", "kfold", "loo"))

    f <- quote(
        function(){

            Atrain <- lapply(bigA, function(x) x[-inds, ])

            if (class(rgcca) == "sgcca")
                tau <- rgcca$c1
            else
                tau <- rgcca$tau
 
            if (rgcca$superblock) {
                Atrain <- Atrain[-length(Atrain)]
                rgcca$C <- NULL
            }

            rgcca_k <- rgcca.analyze(
                Atrain,
                rgcca$C,
                superblock = rgcca$superblock,
                tau = tau,
                ncomp = rgcca$ncomp,
                scheme = rgcca$scheme,
                scale = FALSE,
                type = class(rgcca),
                verbose = FALSE,
                init = rgcca$init,
                bias = rgcca$bias,
                tol = rgcca$tol
            )

            rgcca_k$a <- check_sign_comp(rgcca, rgcca_k$a)

            predict.gcca(
                rgcca_k,
                newA = lapply(bigA, function(x) x[inds, ]),
                type = type,
                fit = fit,
                bloc_to_pred = bloc_to_pred,
                bigA = bigA,
                new_scaled = TRUE
            )
        }
    )

    bigA <- rgcca$blocks
    
    if (validation == "loo")
        v_inds <- 1:nrow(rgcca$blocks[[1]])
    if (validation == "kfold") {
        v_inds <- list()
        for (i in 1:rep) {
            inds <- sample(nrow(rgcca$blocks[[1]]))
            inds <- split(inds, sort(inds %% k))
            v_inds <- c(v_inds, inds)
        }
    }

    if (validation == "test") {
        inds <- sample(nrow(rgcca$blocks[[1]]), size = nrow(rgcca$blocks[[1]]) * 0.3)
        scores <- list(eval(f)())
        preds <- scores$res
    }else{
        scores <- parallel::mclapply(
            v_inds, 
            function(i){
                inds <- i
                eval(f)()
            }, mc.cores = nb_cores
        )
    }
    

    if (validation %in% c("loo", "kfold")) {

        preds <- lapply(
            1:length(rgcca$blocks),
            function(x) Reduce(
                rbind, 
                lapply(
                    scores,
                    function(y) y$pred[[x]]
                    )
                )
            )

        if (validation == "kfold") {
            preds <- lapply(
                1:length(rgcca$blocks),
               function(x) Reduce(
                    rbind,
                    lapply(
                        1:nrow(rgcca$blocks[[1]]),
                        function(y) apply( preds[[x]][ row.names(preds[[x]]) == rownames(rgcca$blocks[[1]])[y], ], 2, mean)
                    )
                )
             )
        }
        
        names(preds) <- names(rgcca$blocks)

    for (x in 1:length(preds))
        row.names(preds[[x]]) <- row.names(rgcca$blocks[[1]])


    }
    
    scores <- mean(unlist(lapply(scores, function(x) x$score)))

    return(list(scores = scores, preds = preds))
}
