#' Creates a table with baseline IDF stats 
#' 
#' This is a summarise operation. for concepts given a set of samples
#' given a dataframe whose grouping defines the concept to calculate
#' IDF stats for
#' 
#' @param groupedDf a dataframe whose grouping defines the "term" for which we calulate the IDF
#' @param sampleVars the column(s) that contains the unique id of a sample, i.e. traditionally a "document" but could be a patient. escaped by vars(...)
#' @param countVar (optional) a field that contains a count. If this is given then it is assumed that the concept & document combinations are unique 
#' @param totalSamples if the data is incomplete (not every document has a concept in it then the expected number fo samples can be specified here)
#' @return a data frame with idf stats for each concept in each group (i.e. document)
#' @import dplyr
#' @export
calculateIdf = function(groupedDf, sampleVars, countVar=NULL, totalSamples=NA) {
	grps = groupedDf %>% groups()
	if (length(grps) == 0) stop("dataframe must be grouped, grouping defines the term for which we are calculating IDF")
	countVar = tryCatch(ensym(countVar),error = function(e) NULL)
	
	# nodes are the unique concepts accross all groups - these are counted to give us:
	# count = number of occurrences in corpus
	# documents = number of documents that each concept appears in (useful for idf)
	
	# f_t - is the number of instances of a given term in the corpus
	# f_td - is the number of instances of a given term in a given document
	# tmp = groupedDf %>% groupwiseCount(sampleVars, countVar, summarise=TRUE) %>% rename(f_t = N, f_td = N_x)
	
	# count here is number of instances of term (grps) in document (sampleVars)
	if (identical(countVar,NULL)) {
		tmp = groupedDf %>% group_by(!!!grps, !!!sampleVars) %>% summarise(count = n())
	} else {
		tmp = groupedDf %>% group_by(!!!grps, !!!sampleVars) %>% summarise(count = sum(!!countVar, na.rm=TRUE))
	}
	
	if (is.na(totalSamples)) {
		totalSamples = groupedDf %>% ungroup() %>% select(!!!sampleVars) %>% distinct() %>% count() %>% pull(n)
	}
	
	nodes = tmp %>% ungroup() %>% group_by(!!!grps) %>% #,isOrigin) %>% 
			# f_t - is the number of instances of a given term in the corpus
			# N_t - is number of documents containing given term in the corpus 
			summarise(
					f_t=sum(count, na.rm=TRUE), 
					N_t=n()
			) %>% ungroup() %>%
			# max_N_t = max number of documents containing any given single term in corpus (term that appears in most documents)
			# max_f_t = max number of any given term occurrences in in corpus (term that appears in most documents)
			# N = total number of documents in the corpus
			mutate(
					max_N_t = max(N_t, na.rm=TRUE),
					max_f_t = max(f_t, na.rm=TRUE),
					N = local(totalSamples)
			) 
	    # %>% compute(unique_index=sapply(grps, as.character))
	
	return(nodes)
}

#' Calculate a term frequecy / inverse document frequency for concepts associated with a grouping
#' 
#' The grouped dataframe here acts as a "Document" from the perpective of the TFIDF calculation but might be a person
#' TODO: fix methods of calculating tfidf
#' 
#' @param groupedDf a dataframe whose grouping defines the "term"
#' @param sampleVars the column(s) that contains the unique id of a sample, i.e. traditionally a "document" but could be a patient. escaped by vars(...). This can include an outcome variable.
#' @param countVar a field that contains a count. If this is given then it is assumed that the concept & document combinations are unique 
#' @param idfDf an optional data frame containing idf information from this or another corpus
#' @param k1 default 1.2 - okapi BM25 parameter
#' @param b default 0.95 - okapi BM25
#' @return a data frame with tfidf stats for each concept in each group (i.e. document)
#' @import dplyr
#' @export
calculateTfidf = function(groupedDf, sampleVars, countVar=NULL, idfDf=NULL, k1 = 1.2, b = 0.95) {
	
	grps = groupedDf %>% groups()
	if (length(grps) == 0) stop("dataframe must be grouped - grouping defines the terms")
	countVar = tryCatch(ensym(countVar),error = function(e) NULL)
	
	termsList = sapply(grps,as.character)
	docsList = sapply(sampleVars,as_label)
	joinList = df %>% joinList(sampleVars)
	
	totalSamples = groupedDf %>% ungroup() %>% select(!!!sampleVars) %>% distinct() %>% count() %>% pull(n)
	if(totalSamples <= 1) {
		stop("Only one document in collection - tfidf is meaningless")
	}
	
	# get the count of codes/terms - n
	# totalCodes = nodes %>% summarise(n=sum(count)) %>% pull(n) %>% collect()
	# f = totalCodes 
	# n.b. this is never used
	
	if (identical(countVar,NULL)) {
		
		if (identical(idfDf,NULL)) {
			idfDf = calculateIdf(groupedDf, sampleVars, totalSamples=totalSamples)
		}
		groupedDf = groupedDf %>% group_by(!!!grps, !!!sampleVars) %>% summarise(count = n())
		
	} else {
		
		if (identical(idfDf,NULL)) {
			idfDf = calculateIdf(groupedDf, sampleVars, !!countVar, totalSamples=totalSamples)
		}
		groupedDf = groupedDf %>% group_by(!!!grps, !!!sampleVars) %>% summarise(count = sum(!!countVar, na.rm=TRUE))
	}
	
	# in groupedDf we are counting concepts / terms
	# the tf statistic
	
	# Doucment level grouping
	# count is number of occurences of given concept(term - t) in each document(d) this is n_td
	# f_td frequency of a given term in a given document (also referred to a f_td)
	# f_d is the total n of codes in a given d
	# max_f_d maximum number of terms in each document
	# max_f_td greatest number of any term in a given document (also referred to a max_f_td)
	groupedDf = groupedDf %>% ungroup() %>% group_by(!!!sampleVars) %>% mutate(
			f_td = count,
			f_d = sum(count, na.rm=TRUE),
			max_f_d=max(f_d, na.rm=TRUE),
			max_f_td = max(count, na.rm=TRUE)
	)
	
	# Corpus level stats
	# avg_n_d average number of terms in all documents - is the same as avg_dl - document length (in terms of codes) on Nr in D
	groupedDf = groupedDf %>% ungroup() %>% mutate(
			avg_f_d=mean(as.double(f_d),na.rm=TRUE),
	)
	
	# G. Paltoglou and M. Thelwall, “A Study of Information Retrieval Weighting Schemes for Sentiment Analysis,” in Proceedings of the 48th Annual Meeting of the Association for Computational Linguistics, Uppsala, Sweden, 2010, pp. 1386–1395 [Online]. Available: http://dl.acm.org/citation.cfm?id=1858681.1858822
	# TODO: Other tfidf metrics
	# TODO: how do we validate this?
	# TODO: split out idf generation so we can do for whole corpus
	
	tmp = groupedDf %>% left_join(
					idfDf %>% select(!!!grps, f_t, N_t, max_N_t, max_f_t, N),
					by=termsList
			) %>% mutate(
					# TODO methods of TFIDF
					length_normalised_tfidf = as.double(f_td) / f_d * log( as.double(N) / N_t),
					o_bm25 = ((k1+1)*as.double(f_td))/(k1*((1-b)+b*as.double(f_d)/avg_f_d)+f_td),
					k_bm25 = log(as.double(N-N_t)+0.5)/(N_t+0.5)
			) %>% mutate (
					okapi_bm25 = o_bm25*k_bm25
			) %>% mutate(
					sum_sq = sqrt(sum(okapi_bm25^2, na.rm=TRUE))
			) %>% mutate(
					norm_okapi_bm25 = okapi_bm25/sum_sq
			)
	
	return(tmp %>% group_by(!!!grps))
}

#' Calculate a co-occurrence matrix for concepts associated with a grouping
#' 
#' The grouped dataframe here acts as a "Document" from the perpective of the co-occurrence but might be a person. For the cose of co-occurrence we are only interested in the existence of
#' a term in a document and not the counts of either.
#' 
#' @param groupedDf a dataframe whose grouping defines the "document" and contains all observations of the terms in that document
#' @param conceptIdVar a field that contains the unique id of a "term"
#' @param filteredIdfDf an optional data frame containing idf information from this or another corpus. It is recommended that this is a heavily filtered list as the result will be the square of the size of this
#' @param bigResult are you expecting a big result?
#' @return a data frame with co-occurrence stats for each concept in filteredIdf (i.e. document)
#' @import dplyr
#' @export
calculateCooccurrenceMI = function(groupedDf, conceptIdVar, filteredIdfDf = NA, bigResult=FALSE) {
	
	grps = groupedDf %>% groups()
	if (length(grps) == 0) stop("dataframe must be grouped")
	conceptIdVar = ensym(conceptIdVar)
	conceptId1Var = as.symbol(paste0(as.character(conceptIdVar),"1"))
	conceptId2Var = as.symbol(paste0(as.character(conceptIdVar),"2"))
	joinList = sapply(grps,as.character)
	#indexList = c(grpsList,as.character(conceptIdVar))
	
	if (identical(filteredIdfDf,NA)) {
	  # calculate the IDF assuming the groupedDf contains the whol corpus
	  filteredIdfDf = calculateIdf(groupedDf, !!conceptIdVar)
	} else {
	  # make sure the IDF does not include concepts outside of the groupedDf set, otherwise these will generate zero co-occurrence entries in the final result
	  filteredIdfDf = filteredIdfDf %>% semi_join(groupedDf, by=as.character(conceptIdVar))
	}
	
	# cut down the groupedDf to just the subset and 
	tmpGroupedDf = groupedDf %>% inner_join(filteredIdfDf, by=as.character(conceptIdVar)) 
	tmpGroupedDf = tmpGroupedDf %>% select(!!!grps, !!conceptIdVar)
	tmpGroupedDf = tmpGroupedDf %>% distinct(!!!grps, !!conceptIdVar) 
	# TODO: for some very obscure reason this was throwing an error when using compute(), complaining about a name collision. The SQL it was generating was fine though, and worked seperatey
	# I think this has eomthing to do with the distinct cause which wraps the sql using a subquery identifier that does not conform to the same syntax as the others.
	filteredGroupedDf = tmpGroupedDf %>% collapse() # %>% compute()
	filteredIdfDf = filteredIdfDf %>% ungroup()
	
	if (filteredIdfDf %>% count() %>% pull(n) > 500 && bigResult == FALSE) stop("The co-occurrence result is going to be greater than 250000. Set bigResult = TRUE to execute.")
	
	# Generate the "full" matrix of cooccurences
	# is actually calculated as MI is symmetrical
	# created on a join using a constant value j to give all possible cooccurrences.
	# add corpus level frequency information from each of the 2 sides of the co-occurrence matrix
	fullMatrix = filteredIdfDf %>% select(tmp_concept_id1 = !!conceptIdVar, N_t1 = N_t) %>% mutate(j=1L) %>% inner_join(
			filteredIdfDf %>% select(tmp_concept_id2 = !!conceptIdVar, N_t2 = N_t, N) %>% mutate(j=1L), by="j"
	) %>% select(-j)
	fullMatrix = fullMatrix %>% compute()
	
	# create the co-occurence information from the document level information
	# create a join specfication based on the "document" groupings as a named vector
	# joins = sapply(grps, as.character)
	# names(joins) = lapply(grps, as.character)
	
	# apply the join on the document level infoto work out the cooccurences
	tmp = filteredGroupedDf %>% rename(tmp_concept_id1 = !!conceptIdVar) %>% 
			inner_join(filteredGroupedDf %>% rename(tmp_concept_id2 = !!conceptIdVar), by=joinList)
	
	# count the number of co-occurences
	tmp = tmp %>% ungroup() %>% group_by(tmp_concept_id1,tmp_concept_id2) %>% summarise(N_t1t2 = n()) %>% 
			compute()
	
	# add document level co-occurrence frequency info to the matrix
	fullMatrix = fullMatrix %>% left_join(tmp, by=c("tmp_concept_id1","tmp_concept_id2")) 
	fullMatrix = fullMatrix %>% mutate(N_t1t2 = ifelse(is.na(N_t1t2),0,N_t1t2)) 
	fullMatrix = fullMatrix %>% rename(!!conceptId1Var := tmp_concept_id1, !!conceptId2Var := tmp_concept_id2 )
	
	return(fullMatrix %>% tidyinfostats::probabilitiesFromCounts(N_t1t2, N_t1, N_t2, N) %>% tidyinfostats::calculateBinaryMI())
	
	
}