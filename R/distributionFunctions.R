#' Distribution class
#' 
#' The distribution class wrapps PDF and Quatile functions from a number of distributions and provides some simple stats
#' for those functions, including a sampling function
#'
#' @keywords distributions
#' @import dplyr
#' @import ggplot2
#' @export
Distribution = R6::R6Class("Distribution", public=list(
				
				#### Fields ----
				#' @field pdf the density function
				pdf = NULL,
				#' @field invCdf the centile function
				invCdf = NULL,
				#' @field dots the parameters of the pdf
				dots = NULL,
				
				#### Methods ----
				#' @description Creates a distribution
				#' @param density a function that accepts at least a vector of x value (e.g. dnorm)
				#' @param quantile a function that accepts at least  a vector of p values (e.g. qnorm)
				#' @param ... passed to pdfFunction and centileFunction
				initialize = function(density, quantile, ...) {
					self$dots = rlang::list2(...)
					self$pdf = density
					self$invCdf = quantile
				},
				
				#' @description calculate the probability density at x
				#' @param x a vector of values of x
				#' @return a vector of values of (p(x))
				p = function(x) {
					rlang::exec(self$pdf,x,!!!self$dots)
				},
				
				#' @description calculate the value of x for a centile y
				#' @param y a vector of centiles
				#' @return a vector of values of x
				q = function(y) {
				  rlang::exec(self$invCdf,y,!!!self$dots)
				},
				
				#' @description the PDF and CDF for a distribution as a dataframe
				#' @param xmin the smallest value
				#' @param xmax the smallest value
				#' @param resolution the number of increments in x
				#' @return a dataframe
				getPdf = function(xmin, xmax, resolution=1001) {
					df = tibble(x=seq(xmin,xmax,length.out = resolution))
					dots = self$dots
					df = df %>% mutate(
							px = self$pdf(x,!!!dots),
							CDFx = cumsum(
									((px+lag(px,default=0))/2)*((xmax-xmin)/resolution)
							)
					) 
					#%>% mutate(
					#  px = px/max(CDFx),
					#  CDFx = CDFx/max(CDFx)
					#)
					return(df)
				},
				
				#' @description defines a inverse CDF as a data frame
				#' @param resolution the number of increments in y
				#' @return a dataframe
				getInverseCdf = function(resolution=1001) {
					df =  tibble(
							y=seq(0,1,length.out = resolution)
					)
					dots = self$dots
					df = df %>% mutate(
							x = self$invCdf(y,!!!dots)
					)
					return(df)
				},
				
				#' @description calculates the integral of -p(x)*log(p(x)) from -Infinity to Infinity
				#' @return a value
				theoreticalEntropy = function() {
				  fn = function(x) ifelse(self$p(x)==0,0,-self$p(x)*log(self$p(x)))
				  H = integrate(fn,-Inf,Inf,subdivisions=2000, rel.tol=.Machine$double.eps^.125)$value
				  return(H)
				},
				
				#' @description gets a label for this distribution based on the parameters passed
				#' @return a string
				label = function() {
					paste(paste0(names(self$dots),"=",twoDp(self$dots)),collapse="; ")
				},
				
				#' @description gets a label for this distribution based on the parameters passed
				#' @param y a vector of y values to apply the quantile function to 
				#' @return a daat frame of y and corresponding x values 
				inverse = function(y) {
					dots = self$dots
					return(data.frame(y=y,x=rlang::exec(self$invCdf,y,!!!dots)))
				},
				
				#' @description produce a set of samples conforming to this distribution
				#' @param n the number of samples
				#' @return a data frame of samples (labelled x)
				sample = function(n=10) {
					dots = self$dots
					return(data.frame(x=rlang::exec(self$invCdf,runif(n),!!!dots)))
				},
				
				#' @description plot this dictributions as pdf and cdf
				#' @param xmin the minimum of the support range to plot
				#' @param xmax the maximum of the support range to plot
				#' @param resolution the number of points to generate for the plot (default 1001)
				#' @return a ggassemble plot object
				plot = function(xmin, xmax, resolution=1001) {
					df = self$getPdf(xmin, xmax, resolution)
					pdfPlot = ggplot(df,aes(x=x,y=px,ymax=px))+geom_line(colour="red")+geom_area(fill="red",alpha=0.3)+ylab("p(x)")
					cdfPlot = ggplot(df,aes(x=x,y=CDFx))+geom_line(colour="blue")+ylab("P(x)")
					return(patchwork::wrap_plots(pdfPlot,cdfPlot,nrow=1))
				},
				
				#' @description a 2 dp formatter
				#' @param x a number
				#' @param unit a unit
				twoDp = function(x,unit="") {
				  paste0(sprintf("%.2f",x),unit)
				}
		
		))

#### -------


#' Get a normal distribution wrapper
#'
#' @keywords distributions
#' @import dplyr
#' @export
NormalDistribution = R6::R6Class("NormalDistribution", inherit=Distribution, public=list(
				#' @field mu the mean of the normal distribuition  
				mu=NULL,
				#' @field sigma the sigma value for the distribution 
				sigma=NULL,
				
				#' @description get a new normal distribution parameterised by mean and sd
				#' @param mean the mean
				#' @param sd the sd
				#' @return a NormalDistribution object
				initialize = function(mean=runif(1,-2,2),sd=runif(1,0.5,5)) {
					self$mu = mean
					self$sigma = sd
					super$initialize(density=dnorm,quantile=qnorm,mean=mean,sd=sd)
				},
				
				#' @description gets a label for this distribution based on the parameters passed
				#' @return a string
				label = function() {
					return(paste0("Norm: \U003BC=",self$twoDp(self$mu),"; \u03C3=",self$twoDp(self$sigma)))
				},
				#' @description calculates the theoretical entropy 
				#' @return a value
				theoreticalEntropy = function() {
				  log(self$sigma*sqrt(2*pi*exp(1)))
				}
		))

#' Get a log normal distribution wrapper
#'
#' @keywords distributions
#' @import dplyr
#' @export
LogNormalDistribution = R6::R6Class("LogNormalDistribution", inherit=Distribution, public=list(
				#' @field mu the mean of the normal distribuition
				mu=NULL,
				#' @field sigma the sd of the normal distribuition
				sigma=NULL,
				#' @description get a LogNormal distribution based on tow paramteretisation options - mean and sd (on natural scale) or mode and sd (on natural scale)
				#' @param mode the mode on a natural scale
				#' @param sd the standard deviation on a natural scale
				#' @param mean the mean on a natural scale
				initialize = function(mode=runif(1,1,4),sd=runif(1,0.5,5),mean=NA) {
					if (is.na(mean)) {
					  self$mu = log( (mode + sqrt(mode^2 + 4 * sd^2))/2  )
				    self$sigma = sqrt(self$mu-log(mode))
					} else {
					  #https://en.wikipedia.org/wiki/Log-normal_distribution LN7 -> LN2
					  self$mu = log(mean/sqrt(1+(sd^2)/(mean^2)))
					  self$sigma = sqrt((log(mean)-self$mu)*2)
					}
					super$initialize(density=dlnorm,quantile=qlnorm,meanlog=self$mu,sdlog=self$sigma)
				},
				#' @description gets a label for this distribution based on the parameters passed
				#' @return a string
				label = function() {
					return(paste0("LogNorm: \U003BC=",self$twoDp(self$mu),"; \u03C3=",self$twoDp(self$sigma)))
				},
				#' @description calculates the theoretical differential entropy
				#' @return a value
				theoreticalEntropy = function() {
				  self$mu+0.5*log(2*pi*exp(1)*self$sigma^2)
				}
		))


#' Get a uniform distribution wrapper
#'
#' @keywords distributions
#' @import dplyr
#' @export
UniformDistribution = R6::R6Class("UniformDistribution", inherit=Distribution, public=list(
        #' @field min the min
        min=NULL,
        #' @field max the max
        max=NULL,
				#' @description get a Uniform distribution wrapper
				#' @param min the min value of the uniform distribution
				#' @param max the max value of the uniform distribution
				initialize = function(min=runif(1,-3,3),max=min+runif(1,0.5,6)) {
				  self$min = min
				  self$max = max
					super$initialize(density=dunif,quantile=qunif,min=min,max=max)
				},
				
				#' @description get a label for this distribution
				label = function() {
				  return(paste0("Unif: ",paste(paste0(names(self$dots),"=",self$twoDp(self$dots)),collapse="; ")))
				},
				#' @description calculates the theoretical entropy of the distribution
				#' @return a value
				theoreticalEntropy = function() {
				  log(self$max-self$min)
				}
		))

#' Get a reversed Kumaraswamy distribution wrapper 
#' 
#' This is a Kumaraswamy distrbution mirrored in the line x=0.5, with support from x=0..1
#'
#' @keywords distributions
#' @import dplyr
#' @export
MirroredKumaraswamyDistribution = R6::R6Class("MirroredKumaraswamyDistribution", inherit=Distribution, public=list(
  #' @field mode the mode of the kumaraswamy
  mode=NULL,
  #' @field iqr the iqr of the normal distribuition
  iqr=NULL,
  #' @description Kumaraswamy distribution
     #' @param mode the mode of the distribution
     #' @param iqr the iqr of the target (must be strictly less than 0.5)
     initialize = function(mode=runif(1,0.1,0.9),iqr=runif(0.1,0.3)) {
       self$mode = mode
       self$iqr = iqr
       fn = function(a,mode,iqr) (1-0.25^((a*mode^a)/(mode^a+a-1)))^(1/a)-(1-0.75^((a*mode^a)/(mode^a+a-1)))^(1/a)-iqr
       a = stats::uniroot(fn, interval=c(1, 10000), mode=(1-mode), iqr=iqr)$root
       b = (-1+a+(1-mode)^a)/(a*(1-mode)^a)
       super$initialize(
         density=function(x,a,b) a*b*(1-x)^(a-1)*(1-(1-x)^a)^(b-1),
         quantile=function(y,a,b) 1-(1-y^(1/b))^(1/a),
         a=a,b=b)
     },
     
     #' @description get a label for this distribution
     label = function() {
       return(paste0("Kum: mode=",self$twoDp(self$mode),"; var=",self$twoDp(self$var)))
     }
))

#' Get a kumaraswamy distribution wrapper
#'
#' This is a Kumaraswamy distrbution with support from x=0..1
#'
#' @keywords distributions
#' @import dplyr
#' @export
KumaraswamyDistribution = R6::R6Class("KumaraswamyDistribution", inherit=Distribution, public=list(
  #' @field mode the mode of the kumaraswamy
     mode=NULL,
     #' @field iqr the iqr of the normal distribuition
     iqr=NULL,
     #' @description Kumaraswamy distribution
     #' @param mode the mode of the distribution
     #' @param iqr the iqr of the target (must be strictly less than 0.5)
     initialize = function(mode=runif(1,0.1,0.9),iqr=runif(0.1,0.3)) {
       self$mode = mode
       self$iqr = iqr
       fn = function(a,mode,iqr) (1-0.25^((a*mode^a)/(mode^a+a-1)))^(1/a)-(1-0.75^((a*mode^a)/(mode^a+a-1)))^(1/a)-iqr
       a = stats::uniroot(fn, interval=c(1, 10000), mode=mode, iqr=iqr)$root
       b = (-1+a+mode^a)/(a*mode^a)
       super$initialize(
         density=function(x,a,b) a*b*x^(a-1)*(1-x^a)^(b-1),
         quantile=function(y,a,b) (1-(1-y)^(1/b))^(1/a),
         a=a,b=b)
     },
     
     #' @description get a label for this distribution
     label = function() {
       return(paste0("Kum: mode=",self$twoDp(self$mode),"; var=",self$twoDp(self$var)))
     }
))


##### ---------------


#' Combine continuous distributions
#' 
#' Get a wrapper for multiple (independent) continuous distributions conditioned on a discrete variable
#'
#' @keywords distributions
#' @import dplyr
#' @export
ConditionalDistribution = R6::R6Class("ConditionalDistribution", public=list(
				
				#### Fields ----
				#' @field classes the class names
				classes = c(),
				#' @field weights the relative weights of each distribution
				weights = c(),
				#' @field dists the distribution as Distribution objects
				dists = c(),
				
				#### Methods ----
				#' @description adds in a Distribution with a class name and weight
				#' @param distribution the pdf as an R6 Distribution object (Distribution$new(fn, fnParams...))
				#' @param class the classname
				#' @param weight the relative weight of this class
				withDistribution = function(distribution,class=distribution$label(),weight=1) {
					self$classes = c(self$classes,class)
					self$weights = c(self$weights,weight)
					self$dists = c(self$dists,distribution)
					invisible(self)
				},
				
				#' @description adds a set of random (uniform, normal, lognormal) distributions to the composite distribution with sensible (random) defaults
				#' @param n the number of distributions to add
				withRandomDistributions = function(n=2) {
					for (i in c(1:n)) {
						switch ( sample(1:3, 1),
								self$withDistribution(UniformDistribution$new(),sample(1:5, 1)),
								self$withDistribution(NormalDistribution$new(),sample(1:5, 1)),
								self$withDistribution(LogNormalDistribution$new(),sample(1:5, 1))
						)
					}
					invisible(self)
				},
				
				#' @description produce a set of samples conforming to this distribution, with random discrete value
				#' @param n the number of samples
				#' @return a data frame of samples (labelled x) associated with classes (labelled "class")
				sample = function(n=1000) {
				  max = sum(self$weights)
				  cutPoints = c(-Inf,cumsum(weights))
				  distIndex = cut(runif(n,0,max),cutPoints,labels=FALSE,right=FALSE)
				  classList = self$classes[distIndex]
				  return(sampleByClass(classList))
				},
				
				
				#' @description produce a set of samples conforming to this distribution with preset discrete value
				#' @param classVector the number of samples
				#' @return a data frame of samples with continuous values (labelled x) associated with discrete classes (labelled "y") and a sample id column (labelled i)
				sampleByClass = function(classVector) {
				  out = NULL
				  id = 1
				  for(class in classVector) {
				    dist = self$dists[self$classes==class][[1]]
				    out = out %>% bind_rows(tibble(
				      i = id,
				      x = dist$sample(1)$x,
				      y = class))
				    id = id+1
				  }
				  return(out)
				},
				
				
				#' @description get the pdf of these distributions as a dataframe
				#' @param xmin - the minimum of the support
				#' @param xmax - the maximum of the support
				#' @param resolution -  the number of points in the pdf
				#' @return a datafram containing all classes
				getPdf = function(xmin, xmax, resolution=1001) {
					out = NULL
					py = self$weights/sum(self$weights)
					for(i in c(1:length(self$classes))) {
						out = out %>% rbind(
								self$dists[[i]]$getPdf(xmin,xmax,resolution) %>% mutate(
										y = self$classes[[i]], 
										py=py[[i]], 
										pxy=px*py[[i]], 
										CDFxy=CDFx*py[[i]]))
					}
					out = out %>% group_by(x) %>% mutate(px = sum(pxy), CDFx = sum(CDFxy))
					return(out)
				},
				
				#' @description plot this distributions as pdf and cdf
				#' @param xmin - the minimum of the support
				#' @param xmax - the maximum of the support
				#' @return a ggassemble plot object
				plot = function(
				    xmin=self$theoreticalMean()-3*sqrt(self$theoreticalVariance()),
				    xmax=self$theoreticalMean()+3*sqrt(self$theoreticalVariance())
				) {
					out = self$getPdf(xmin=xmin,xmax=xmax)
					pdfPlot = ggplot(out,aes(x=x,y=pxy,ymax=CDFxy,colour=y,fill=y))+geom_line(aes(y=px),colour="grey50")+geom_line()+geom_area(alpha=0.3,position="identity")+ylab("p(x\u2229y)")+guides(fill="none", colour="none")#+theme(legend.position = "bottom",legend.direction = "vertical")
					cdfPlot = ggplot(out,aes(x=x,y=CDFxy,colour=y))+geom_line(aes(y=CDFx),colour="grey50")+geom_line()+ylab("P(x\u2229y)")#+theme(legend.position = "bottom",legend.direction = "vertical")
					return(patchwork::wrap_plots(pdfPlot,cdfPlot,nrow=1,guides="collect"))
				},
				
				#' @description generate the theoretical mutual information for this set of distributions using numerical integration of the underlying functions
				#' @return a single value for the mutual information of this function
				theoreticalMI = function() {
					py = self$weights/sum(self$weights) # y is classes
					py_matrix = function(x) matrix(rep(py,length(x)),nrow = length(x),byrow=TRUE)
					pxy = function(x) py_matrix(x)*sapply(self$dists, function(dist) dist$p(x)) # i.e. p(y) * p(x given y)
					px = function(x) rowSums(pxy(x)) # this gives us over all x
					px_matrix = function(x) matrix(rep(px(x),length(py)),ncol=length(py))
					
					Iix = function(x) rowSums(pxy(x)*log( pxy(x) / (px_matrix(x)*py_matrix(x))),na.rm = TRUE)
					# Iix = function(x) sapply(x, function(x) sum( pxy(x)*log(pxy(x)/(py*px(x)) ), na.rm = TRUE))
					# browser()
					#TODO: can throw error here due to non finite function. To investigate.
					return(integrate(Iix,-Inf,Inf,subdivisions = 2000, rel.tol=.Machine$double.eps^.125)$value)
				},
				
				#' @description generate the theoretical mean
				#' @return a single value for the mean of this function
				theoreticalMean = function() {
					py = self$weights/sum(self$weights) # y is classes
					pxy = function(x) py*sapply(self$dists, function(dist) dist$p(x)) # i.e. p(y) * p(x given y)
					px = function(x) sapply(x, function(x) sum(pxy(x))) # this gives us over all x
					xpx = function(x) x*px(x)
					return(integrate(xpx,-Inf,Inf, subdivisions=10000)$value)
				},
				
				#' @description generate the theoretical variance
				#' @return a single value for the variance of this function
				theoreticalVariance = function() {
					mu = self$theoreticalMean()
					py = self$weights/sum(self$weights) # y is classes
					pxy = function(x) py*sapply(self$dists, function(dist) dist$p(x)) # i.e. p(y) * p(x given y)
					px = function(x) sapply(x, function(x) sum(pxy(x))) # this gives us over all x
					sdpx = function(x) ((x-mu)^2)*px(x)
					return(integrate(sdpx,-Inf,Inf, subdivisions=10000)$value)
				}
		
		))

##### ---------------


#' Simulate a data set with multiple features
#' 
#' simulating a dataset with more than one feature requires some logic sampling
#'
#' @keywords distributions
#' @import dplyr
#' @export
MultivariableDistribution = R6::R6Class("MultivariableDistribution", public=list(
  
  #### Fields ----
  #' @field features the feature names
  features = c(),
  #' @field classes the class names
  classes = c(),
  #' @field weights the relative weights of each distribution
  weights = c(),
  #' @field dists the distribution as Distribution objects
  dists = c(),
  
  #### Methods ----
  #' @description adds in a ConditionalDistribution with a class name
  #' @param featureName the class name
  #' @param distribution the pdf as an R6 ConditionalDistribution object (ConditionalDistribution$new(fn, fnParams...))
  withConditionalDistribution = function(distribution,featureName) {
    self$features = c(self$features,featureName)
    self$dists = c(self$dists,distribution)
    invisible(self)
  },
  
  #' @description sets the relative weights of the different outcome classes in the simulation
  #' @param listWeights the weights as a named list e.g. list(feature1 = 0.1, ...)
  withClassWeights = function(listWeights) {
    self$classes = names(listWeights)
    self$weights = unlist(listWeights, use.names = FALSE)
    sapply(self$dists, function(cls) cls$weights = unlist(listWeights[cls$classes]))
    invisible(self)
  },
  
  #' @description produce a set of samples conforming to these distributions
  #' @param n the number of samples
  #' @return a data frame of samples (labelled x) associated with classes (labelled "class")
  sample = function(n=1000) {
    # TODO: change this to allow for missing values to be produced based on a class by class basis decided by the weights of the conditional disctributions
    max = sum(self$weights)
    cutPoints = c(-Inf,cumsum(self$weights))
    distIndex = cut(runif(n,0,max),cutPoints,labels=FALSE,right=FALSE)
    classList = self$classes[distIndex]
    out = NULL
    for (i in c(1:length(self$dists))) {
      dist = self$dists[[i]]
      # dist$weights
      # runif(n,0,1) < 
      featureName = self$features[[i]]
      tmp = dist$sampleByClass(classList)
      tmp = tmp %>% mutate(feature=featureName)
      out = out %>% bind_rows(tmp)
    }
    return(out %>% rename(outcome=y,value=x,sample=i))
  },
  
  #' @description plot this distributions as pdf and cdf
  #' @param xmin - the minimum of the support
  #' @param xmax - the maximum of the support
  #' @return a ggassemble plot object
  plot = function() {
    plots = lapply(self$dists, function(d) d$plot())
    return(patchwork::wrap_plots(plots,ncol=1,guides="collect"))
  }
  
))