img.Test <- function(batch,pset,x) {
par(mfrow = c(2,2))
image(batch[,x])
image(pset, type = "weights", which = x)
image(pset, type = "resids", which = x)
image(pset, type = "sign.resids", which = x)
}