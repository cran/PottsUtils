SW <- function(n, nvertex, ncolor, edges, weights=1, beta){
    if (ncol(edges) != 2)
        stop("'edges' with two columns have to be provided.")
    nedge <- nrow(edges)
    if (nedge < length(weights))
        stop("The number of 'edges' is less than the number of 'weights'.")
    if (length(weights) < nedge)
        weights <- rep(weights, length=nrow(edges))

    oneIteration <- sample(x=1:ncolor, nvertex, replace=TRUE)
    colors <- matrix(0, nrow=nvertex, ncol=n)
    for(i in 1:n){
        bondsTest <- weights *
            ifelse(oneIteration[edges[,1]] - oneIteration[edges[,2]] == 0, 1, 0)
        bondsProb <- 1 - exp(-beta * bondsTest)
        bonds <- edges[runif(nedge) < bondsProb,]
        patches <- getPatches(bonds, nvertex)
        npatch <- length(patches)
        newColors <- sample(x=1:ncolor, npatch, replace=TRUE)
        for (j in 1: npatch){
            oneIteration[patches[[j]]] <- newColors[j]
        }
        colors[,i] <- oneIteration
    }
    colors
}

