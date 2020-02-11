#' ClassifierResult class
#' 
#' result
#'
#' @keywords distributions
#' @import dplyr
#' @import ggplot2
#' @export
ClassifierResult = R6::R6Class("ClassifierResult", public=list(
  
  #### Fields ----
  #' @field data - the classifier prediction probabilities and outcomes
  data = NULL,
  #' @field classes - the classes of the outcome as a vector
  classes = NULL,

  #### Methods ----
  #' @description Creates a distribution
  #' @param predictionProbabilities a matrix of probabilities with rows for every sample, and colums for every possible class
  #' @param actual the gold standard output as a vector of outcomes for each sample
  #' @param ... passed to pdfFunction and centileFunction
  initialize = function(data, classes) {
    self$data = data
    self$classes = classes
  },
  
  getInfoStats = function() {
    tmp = self$data %>% group_by(prediction) %>% arrange(probability) %>% 
      mutate(
        rank = row_number() # predicted negative (for a given cut off)
      ) %>% 
      mutate(
        N=max(rank), # total
        match = ifelse(prediction==obs,1,0)
      ) %>% mutate(
        matches=cumsum(match), # false negatives
        total_matches = sum(match) # total actual positives
      ) %>% mutate(
        # at a given value of probability we have
        truePos = total_matches-matches,
        predPos = N-rank, # total - predicted negative
        obsPos = total_matches,
        total = N
      ) %>% select(
        -c(matches,rank,N,total_matches)
      )
    tmp = tmp %>% probabilitiesFromCounts(truePos,obsPos,predPos,total) %>% calculateConfusionMatrixStats()
    return(tmp)
  },
  
  plotRoc = function() {
    tmp = self$getInfoStats()
    return(
      ggplot(tmp, aes(x=false_pos_rate, y=true_pos_rate, colour=prediction))+
        geom_path()+
        xlab("False positive rate")+
        ylab("True positive rate")+
        geom_segment(x=0,y=0, xend=1, yend=1, colour="grey75",linetype="dashed")+
        expand_limits(x=c(0,1),y=c(0,1))
    )
  }
  
))


#### Multiclass classifier ----

#' MulticlassClassifierResult class
#' 
#' result
#'
#' @keywords distributions
#' @import dplyr
#' @import ggplot2
#' @export
MulticlassClassifierResult = R6::R6Class("MulticlassClassifierResult", inherit=ClassifierResult, public=list(
  
  #### Methods ----
  #' @description Creates a distribution
  #' @param ... passed to pdfFunction and centileFunction
  initialize = function(data, classes) {
    super$initialize(data,classes)
  },
  
  getBinaryClassifiers = function() {
    out = lapply(self$classes, function(c) {
      BinaryClassifierResult$new(
        self$data %>% filter(prediction == c) %>% mutate(obs = ifelse(obs==c,obs,"other")), c
      )
    })
    names(out) <- self$classes
    return(out)
  },
  
  getClasses = function() {
    return(self$classes)
  }
  
))

#### Binary classifier ----

#' ClassifierResult class
#' 
#' result
#'
#' @keywords distributions
#' @import dplyr
#' @import ggplot2
#' @export
BinaryClassifierResult = R6::R6Class("BinaryClassifierResult", inherit=ClassifierResult, public=list(
  
  #### Methods ----
  #' @description Creates a distribution
  #' @param predictionProbabilities a matrix of probabilities with rows for every sample, and colums for every possible class
  #' @param actual the gold standard output as a vector of outcomes for each sample
  #' @param ... passed to pdfFunction and centileFunction
  initialize = function(data, classes) {
    super$initialize(data,classes)
  }
  
))

#' ClassifierResult factory
#' 
#'
#' @keywords distributions
#' @param predictionProbabilities a matrix of probabilities with rows for every sample, and colums for every possible class
#' @param obs the gold standard output as a vector of outcomes for each sample
#' @export
ClassifierResult$fromPredictions = function(predictionProbabilities, obs) {
    tmp = data.frame(predictionProbabilities) %>% mutate(tmp_obs=as.character(obs), sample = row_number())
    tmp = tmp %>% tidyr::pivot_longer(cols = colnames(predictionProbabilities), names_to = "prediction", values_to = "probability") %>% rename(obs = tmp_obs)
    classes = unique(c(tmp$prediction, tmp$obs))
    if(length(classes) < 2) stop("not enough classes...?")
    if(length(classes) == 2) {
      tmp = tmp %>% filter(prediction=classes[[1]])
      return(BinaryClassifierResult$new(tmp,classes))
    } else {
      return(MulticlassClassifierResult$new(tmp,classes))
    }
}
