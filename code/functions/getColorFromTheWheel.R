# getColorFromTheWheel ===========
# select colors from the color wheel for a given
#=================================
getColorFromTheWheel <- function(f.scores, design, order.dim = 1, color.range = c("#006e00", "#efc900")){
  mean.f <- getMeans(f.scores, design)[,1:2]
  ordered.f <- mean.f[order(mean.f[,1]),]
  colfunc <- colorRampPalette(color.range)
  colvec <- colfunc(nrow(mean.f))
  names(colvec) <- rownames(ordered.f)
  return(colvec)
}