#' Get a distribution wrapper
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
				#' @description Sets up an omop database connection
				#' @param density a function that accepts at least an x value (e.g. dnorm)
				#' @param quantile a function that accepts at least an y value (e.g. qnorm)
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
				
				#' @description calculates the integral of -p(x)*log(p(x)) from -Infinity to infinity
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
				#' @field sigma the  
				sigma=NULL,
				
				#' @description plot this dictributions as pdf and cdf
				#' @param mean the mean
				#' @param sd the sd
				#' @return a ggassemble plot object
				initialize = function(mean=runif(1,-2,2),sd=runif(1,0.5,5)) {
					self$mu = mean
					self$sigma = sd
					super$initialize(density=dnorm,quantile=qnorm,mean=mean,sd=sd)
				},
				
				#' @description gets a label for this distribution based on the parameters passed
				#' @return a string
				label = function() {
					return(paste0("Norm: \U003BC=",twoDp(self$mu),"; \u03C3=",twoDp(self$sigma)))
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
				#' @description plot this dictributions as pdf and cdf
				#' @param mode the mode - must be 1 or greater TODO: find out why
				#' @param sd the sd
				initialize = function(mode=runif(1,1,4),sd=runif(1,0.5,5)) {
					self$mu = mode
					self$sigma = sd
					fn = function(mean,mode,sd) exp(mean)*mode-exp(2*mean)+sd^2
					meanlog = stats::uniroot(fn, interval=c(-100,100), mode=mode, sd=sd)$root
					sdlog = sqrt(meanlog-log(mode))
					super$initialize(density=dlnorm,quantile=qlnorm,meanlog=meanlog,sdlog=sdlog)
				},
				#' @description gets a label for this distribution based on the parameters passed
				#' @return a string
				label = function() {
					return(paste0("LogNorm: mode=",twoDp(self$mu),"; var=",twoDp(self$sigma)))
				}
		))


#' Get a uniform distribution wrapper
#'
#' @keywords distributions
#' @import dplyr
#' @export
UniformDistribution = R6::R6Class("UniformDistribution", inherit=Distribution, public=list(
				#' @description plot this dictributions as pdf and cdf
				#' @param min the min value of the uniform distribution
				#' @param max the max value of the uniform distribution
				initialize = function(min=runif(1,-3,3),max=min+runif(1,0.5,6)) {
					super$initialize(density=dunif,quantile=qunif,min=min,max=max)
				},
				
				#' @description get a label for this distribution
				label = function() {
				  return(paste0("Unif: ",paste(paste0(names(self$dots),"=",twoDp(self$dots)),collapse="; ")))
				}
		))

twoDp = function(x,unit="") {
	paste0(sprintf("%.2f",x),unit)
}

##### ---------------


#' Get a named class distribution wrapper
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
				#' @param class the classname
				#' @param weight the relative weight of this class
				#' @param distribution the pdf as an R6 Distribution object (Distribution$new(fn, fnParams...))
				withDistribution = function(weight,distribution,class=distribution$label()) {
					self$classes = c(self$classes,class)
					self$weights = c(self$weights,weight)
					self$dists = c(self$dists,distribution)
					invisible(self)
				},
				
				#' @description adds a set of random distributions to the composite 
				#' @param n the number of distributions to add
				withRandomDistributions = function(n=2) {
					for (i in c(1:n)) {
						switch ( sample(1:3, 1),
								self$withDistribution(sample(1:5, 1), UniformDistribution$new()),
								self$withDistribution(sample(1:5, 1), NormalDistribution$new()),
								self$withDistribution(sample(1:5, 1), LogNormalDistribution$new())
						)
					}
					invisible(self)
				},
				
				#' @description produce a set of samples conforming to this distribution
				#' @param n the number of samples
				#' @return a data frame of samples (labelled x) associated with classes (labelled "class")
				sample = function(n=1000) {
					counts = floor(self$weights*n/sum(self$weights))
					while (sum(counts) < n) {
					  i = sample(length(self$weights),size=1)
					  counts[i] = counts[i]+1
					}
					out = NULL
					for(i in c(1:length(self$classes))) {
						out = out %>% rbind(self$dists[[i]]$sample(counts[[i]]) %>% mutate(y = self$classes[[i]]))
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
				plot = function(xmin,xmax) {
					out = self$getPdf(xmin=xmin,xmax=xmax)
					pdfPlot = ggplot(out,aes(x=x,y=pxy,ymax=CDFxy,colour=y,fill=y))+geom_line(aes(y=px),colour="grey50")+geom_line()+geom_area(alpha=0.3,position="identity")+ylab("p(x\u2229y)")+theme(legend.position = "bottom",legend.direction = "vertical")
					cdfPlot = ggplot(out,aes(x=x,y=CDFxy,colour=y))+geom_line(aes(y=CDFx),colour="grey50")+geom_line()+ylab("P(x\u2229y)")+theme(legend.position = "bottom",legend.direction = "vertical")
					return(patchwork::wrap_plots(pdfPlot,cdfPlot,nrow=1))
				},
				
				# testing the MI generation with a dataframe rather than a function
#				numericalMI = function(xmin,xmax,resolution=1001) {
#					d = self$getPdf(xmin,xmax,resolution)
#					d = d %>% mutate(
#									Ixy = pxy*log(pxy/(px*py)),
#							) %>% group_by(y) %>% arrange(y,x) %>% mutate(
#									# integrate x for each y
#									Idxy = (x-lag(x,default=xmin))*(Ixy+lag(Ixy,default=0))/2
#							) 
#					d2 = d %>% summarise(
#							Iy = sum(Idxy,na.rm=TRUE)
#					)
#					I = d2 %>% summarise(
#							# sum over y
#							I = sum(Iy,na.rm=TRUE)
#					) %>% pull(I)
#					return(I)
#				},
				
				#' @description generate the theoretical mutual information
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
					return(integrate(xpx,-Inf,Inf)$value)
				},
				
				#' @description generate the theoretical variance
				#' @return a single value for the variance of this function
				theoreticalVariance = function() {
					mu = self$theoreticalMean()
					py = self$weights/sum(self$weights) # y is classes
					pxy = function(x) py*sapply(self$dists, function(dist) dist$p(x)) # i.e. p(y) * p(x given y)
					px = function(x) sapply(x, function(x) sum(pxy(x))) # this gives us over all x
					sdpx = function(x) ((x-mu)^2)*px(x)
					return(integrate(sdpx,-Inf,Inf)$value)
				}
		
		))