
#' calculate mutual information from counts
#'
#' @param x1y1 cooccurrence of x1 and y1 in corpus.
#' @param x1 occurrence of x1 in corpus
#' @param y1 occurrence of y1 in corpus
#' @param total total size of corpus
#' @return the mutual information
#' @export  
miFromCounts = function(x1y1, x1, y1, total) {
  notValid = (x1y1 > x1 | x1y1 > y1 | x1 > total | y1 > total);
  x1y0 = x1-x1y1;
  x0y1 = y1-x1y1;
  p_x1y1 = (as.double(x1y1))/total;
  p_x0y1 = (as.double(x0y1))/total;
  p_x1y0 = (as.double(x1y0))/total;
  return(ifelse(notValid,NaN,
                mi(p_x1y1,p_x0y1,p_x1y0)));
}

#' A calculate pointwise mutual information from counts
#'
#' @param x1y1 cooccurrence of x1 and y1 in corpus.
#' @param x1 occurrence of x1 in corpus
#' @param y1 occurrence of y1 in corpus
#' @param total total size of corpus
#' @return the pointwise mutual information
#' @export
pmiFromCounts = function(x1y1, x1, y1, total) {
  notValid = (x1y1 > x1 | x1y1 > y1 | x1 > total | y1 > total);
  return(ifelse(notValid,NA,pmi(
    (as.double(x1y1))/total,
    (as.double(x1))/total,
    (as.double(y1))/total)));
}

#' A calculate pointwise mutual information from counts
#'
#' @param x1y1 cooccurrence of x1 and y1 in corpus.
#' @param x1 occurrence of x1 in corpus
#' @param y1 occurrence of y1 in corpus
#' @param total total size of corpus
#' @return the normalised pointwise mutual information
#' @export  
npmiFromCounts = function(x1y1, x1, y1, total) {
  notValid = (x1y1 > x1 | x1y1 > y1 | x1 > total | y1 > total);
  return(ifelse(notValid,NA,npmi(
    (as.double(x1y1))/total,
    (as.double(x1))/total,
    (as.double(y1))/total)));
}

#' A calculate pointwise mutual information from probabilities
#'
#' @param p_xy probability of cooccurrence of x and y.
#' @param p_x probability of occurrence of x
#' @param p_y probability of occurrence of y
#' @return the pointwise mutual information
#' @export 
pmi = function(p_xy, p_x, p_y) {
  notValid = (p_xy < 0.0 | p_xy > 1.0 | p_x < 0.0 | p_x > 1.0 | p_y < 0.0 | p_y > 1.0 | p_xy > p_x + 0.00001 | (p_xy > p_y + 0.00001 ) );
  returnZero = (p_x==0 | p_y==0); 
  returnNegInf = (p_xy == 0 & p_x > 0.0 & p_y > 0.0);
  return(
    ifelse(notValid, NaN,
           ifelse(returnZero,0,
                  ifelse(returnNegInf, -Inf,
                         log(p_xy/(p_x*p_y))
                  ))));
}

#' A calculate normalised pointwise mutual information from probabilities
#'
#' @param p_xy probability of cooccurrence of x and y.
#' @param p_x probability of occurrence of x
#' @param p_y probability of occurrence of y
#' @return the normalised pointwise mutual information
#' @export 
npmi = function(p_xy, p_x, p_y) {
  notValid = (p_xy < 0.0 | p_xy > 1.0 | p_x < 0.0 | p_x > 1.0 | p_y < 0.0 | p_y > 1.0 | p_xy > p_x + 0.00001 | p_xy > p_y + 0.00001 )
  returnNeg1 = (p_xy == 0.0 & p_x > 0.0 & p_y > 0.0)
  return(
    ifelse(notValid,NaN,
           ifelse(returnNeg1,-1,
                  pmi(p_xy,p_x,p_y)/(-log(p_xy))
           )));
}



#' calculate mutual information from probabilities
#'
#' @param p_x1y1 probability of cooccurrence of x and y.
#' @param p_x0y1 probability of occurrence of y without x
#' @param p_x1y0 probability of occurrence of x without y
#' @return the mutual information
#' @export 
mi = function(p_x1y1, p_x0y1, p_x1y0) {
  notValid = (p_x1y1 < 0.0 | p_x1y1 > 1.0 | p_x0y1 < 0.0 | p_x0y1 > 1.0 | p_x1y0 < 0.0 | p_x1y0 > 1.0 );
  p_x0y0 = 1.0-(p_x1y1+p_x0y1+p_x1y0);
  p_x1 = p_x1y0+p_x1y1;
  p_x0 = 1.0-p_x1;
  p_y1 = p_x1y1+p_x0y1;
  p_y0 = 1.0-p_y1;
  return(ifelse(notValid,NaN,
                (
                  ifelse(p_x1y1==0,0,p_x1y1*pmi(p_x1y1,p_x1,p_y1))+
                    ifelse(p_x0y1==0,0,p_x0y1*pmi(p_x0y1,p_x0,p_y1))+
                    ifelse(p_x1y0==0,0,p_x1y0*pmi(p_x1y0,p_x1,p_y0))+
                    ifelse(p_x0y0==0,0,p_x0y0*pmi(p_x0y0,p_x0,p_y0))
                )
  ))
}

# TODO:
# accuracy, sensitivity, specificity, etc...
