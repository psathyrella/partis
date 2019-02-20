print.fit.bd <- function (x, ...)
{
   nn <- names(x)
   ll <- length(x)
   if (length(nn) != ll)
       nn <- paste("Component", seq.int(ll))
   for (i in seq_len(ll)) {
       if ((nn[i] != "f.lamb") && (nn[i] != "f.mu"))
       {
       cat(nn[i], ":\n")
       print(x[[i]], ...)
       cat("\n")
        }
        else{
            }
   }
   invisible(x)
}