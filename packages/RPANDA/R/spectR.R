spectR<-function (phylo, method = c("standard")) 
{
    skewness <- function(x, na.rm = FALSE) {
        if (is.matrix(x)) 
            apply(x, 2, skewness, na.rm = na.rm)
        else if (is.vector(x)) {
            if (na.rm) 
                x <- x[!is.na(x)]
            n <- length(x)
            (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
        }
        else if (is.data.frame(x)) 
            sapply(x, skewness, na.rm = na.rm)
        else skewness(as.vector(x), na.rm = na.rm)
    }
    sigma = 0.1
    gKernel <- function(x) 1/(sigma * sqrt(2 * pi)) * exp(-(x^2)/2 * 
        sigma^2)
    kernelG <- function(x, mean = 0, sd = 1) dnorm(x, mean = mean, 
        sd = sd)
    dens <- function(x, bw = bw.nrd0, kernel = kernelG, n = 4096, 
        from = min(x) - 3 * sd, to = max(x) + 3 * sd, adjust = 1, 
        ...) {
        if (has.na <- any(is.na(x))) {
            x <- na.omit(x)
            if (length(x) == 0) 
                stop("no finite or non-missing data!")
        }
        sd <- (if (is.numeric(bw)) 
            bw[1]
        else bw(x)) * adjust
        X <- seq(from, to, len = n)
        M <- outer(X, x, kernel, sd = sd, ...)
        structure(list(x = X, y = rowMeans(M), bw = sd, call = match.call(), 
            n = length(x), data.name = deparse(substitute(x)), 
            has.na = has.na), class = "density")
    }
    integr <- function(x, f) {
        if (!is.numeric(x)) {
            stop("The variable of integration \"x\" is not numeric.")
        }
        if (!is.numeric(f)) {
            stop("The integrand \"f\" is not numeric.")
        }
        if (length(x) != length(f)) {
            stop("The lengths of the variable of integration and the integrand do not match.")
        }
        n = length(x)
        integral = 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] + 
            f[1:(n - 1)]))
        return(integral)
    }
    if (method == "standard") {
        e = eigen(graph.laplacian(graph.adjacency(data.matrix(dist.nodes(phylo)), 
            weighted = T), normalized = F), symmetric = T, only.values = T)
        x = subset(e$values, e$values >= 1e-10)
        d = dens(log(x))
        dsc = d$y/(integr(d$x, d$y))
        principal_eigenvalue <- max(x)
        skewness <- skewness(dsc)
        peak_height <- max(dsc)
        gaps <- abs(diff(x))
        gapMat <- as.matrix(gaps)
        modalities <- c(1:length(gapMat))
        gapMatCol <- cbind(modalities, gapMat)
        eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 
            2]))
        res <- list(eigenvalues = x, principal_eigenvalue = principal_eigenvalue, 
            asymmetry = skewness, peakedness = peak_height, eigengap = eigenGap[, 
                1])
    }
    if (method == "normal") {
        e = eigen(graph.laplacian(graph.adjacency(data.matrix(dist.nodes(phylo)), 
            weighted = T), normalized = T), symmetric = T, only.values = T)
        x = subset(e$values, e$values >= 0)
        d = dens(log(x))
        dsc = d$y/(integr(d$x, d$y))
        principal_eigenvalue <- max(x)
        skewness <- skewness(dsc)
        peak_height <- max(dsc)
        gaps <- abs(diff(x))
        gapMat <- as.matrix(gaps)
        modalities <- c(1:length(gapMat))
        gapMatCol <- cbind(modalities, gapMat)
        eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 
            2]))
        res <- list(eigenvalues = x, principal_eigenvalue = principal_eigenvalue, 
            asymmetry = skewness, peakedness = peak_height, eigengap = eigenGap[, 
                1])
    }
    class(res) <- "spectR"
    return(res)
}
