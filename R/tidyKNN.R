#' calculates an approximate KNN for sourceDf entries in targetDf entries
#' 
#' relise on RANN library for actual calculation. Simply reformats the results into a mapping table format
#' 
#' @param sourceDf a dataframe with a unique Id column
#' @param targetDf a dataframe with a unique Id column
#' @param sourceIdVar the unique Id column of the sourceDf
#' @param targetIdVar the unique Id column of the targetDf
#' @param k the k of knn
#' @param matchVars the columns to calculate the distance on - must be numerics and present in both sourceDf and targetDf, escaped by vars(...)
#' @param distanceVar what will the distance measure be named?
#' @param ... other parameters passed onto RANN::nn2()
#' @return a dataframe containing sourceIdVar, targetIdVar, distanceVar and k columns
#' @export
findKNN = function(sourceDf, targetDf, sourceIdVar, targetIdVar, k=4, matchVars=vars(long,lat), distanceVar="distance", rankVar="rank", ...) {
  sourceIdVar = ensym(sourceIdVar)
  targetIdVar = ensym(targetIdVar)
  distanceVar = ensym(distanceVar)
  rankVar = ensym(rankVar)
  if (as.character(sourceIdVar)==as.character(targetIdVar)) stop("sourceIdVar and targetIdVar must have different names")
  
  sourceDf = sourceDf %>% collect() %>% rename(sourceId = !!sourceIdVar)
  targetDf = targetDf %>% collect() %>% rename(targetId = !!targetIdVar)
  
  matches = RANN::nn2(targetDf %>% select(!!!matchVars), sourceDf %>% select(!!!matchVars), k, ...)
  
  # reformat target identifiers
  tmpIdMatch = matches$nn.idx
  colnames(tmpIdMatch) = c(1:k)
  tmpIdMatch = as.data.frame(tmpIdMatch)
  tmpIdMatch = tmpIdMatch %>% mutate(sourceId = sourceDf$sourceId)
  tmpIdMatch = tmpIdMatch %>% pivot_longer(cols = c(1:k), names_to = "k", values_to = "tmpTargetId") %>% mutate(k=as.integer(k)) 
  tmpIdMatch = tmpIdMatch %>% mutate(!!targetIdVar := targetDf$targetId[tmpTargetId]) %>% select(-tmpTargetId)
  
  # reformat distances
  tmpDist = matches$nn.dists
  colnames(tmpDist) = c(1:k)
  tmpDist = as.data.frame(tmpDist)
  tmpDist = tmpDist %>% mutate(sourceId = sourceDf$sourceId)
  tmpDist = tmpDist %>% pivot_longer(cols = c(1:k), names_to = "k", values_to = as.character(distanceVar)) %>% mutate(k=as.integer(k))
  
  idMatch = tmpIdMatch %>% left_join(tmpDist, by=c("sourceId","k")) 
  idMatch = idMatch %>% rename(!!sourceIdVar := sourceId, !!rankVar := k)
  return(idMatch)
}