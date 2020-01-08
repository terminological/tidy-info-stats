#' Calculates multiple 2 class mutual information scores from confusion matrix probabilities in dplyr friendly manner
#'
#' The purpose of this is to make it possible to calculate MI in a DBPLYR sql table
#' 
#' @param df a dataframe containing one observation per row & p_x1y1, p_x0y1, p_x1y0, and p_x0y0 columns (see probabilitiesFromCounts)
#' @return the datatable with additional columns for MI; PMI0 and PMI1; for all various combinations of outcome
#' @import dplyr
#' @export
calculateBinaryMI = function(df) {
	return(
		df %>% mutate(
			pmi_x1y1 = ifelse( p_x1y1==0, ifelse(p_x1==0 | p_y1==0, 0, NA), log(p_x1y1/(p_x1*p_y1)) ),
			pmi_x0y1 = ifelse( p_x0y1==0, ifelse(p_x0==0 | p_y1==0, 0, NA), log(p_x0y1/(p_x0*p_y1)) ),
			pmi_x1y0 = ifelse( p_x1y0==0, ifelse(p_x1==0 | p_y0==0, 0, NA), log(p_x1y0/(p_x1*p_y0)) ),
			pmi_x0y0 = ifelse( p_x0y0==0, ifelse(p_x0==0 | p_y0==0, 0, NA), log(p_x0y0/(p_x0*p_y0)) )
		) %>% mutate(
			I_xy = (
				ifelse(p_x1y1==0|p_x1==0|p_y1==0, 0, p_x1y1*pmi_x1y1)+
				ifelse(p_x0y1==0|p_x0==0|p_y1==0, 0, p_x0y1*pmi_x0y1)+
				ifelse(p_x1y0==0|p_x1==0|p_y0==0, 0, p_x1y0*pmi_x1y0)+
				ifelse(p_x0y0==0|p_x0==0|p_y0==0, 0, p_x0y0*pmi_x0y0)
			)
		) %>% mutate(
			npmi_x1y1 = ifelse( p_x1y1==0, ifelse(p_x1==0 | p_y1==0, 0, -1), pmi_x1y1 / (-log(p_x1y1)) ),
			npmi_x0y1 = ifelse( p_x0y1==0, ifelse(p_x0==0 | p_y1==0, 0, -1), pmi_x0y1 / (-log(p_x0y1)) ),
			npmi_x1y0 = ifelse( p_x1y0==0, ifelse(p_x1==0 | p_y0==0, 0, -1), pmi_x1y0 / (-log(p_x1y0)) ),
			npmi_x0y0 = ifelse( p_x0y0==0, ifelse(p_x0==0 | p_y0==0, 0, -1), pmi_x0y0 / (-log(p_x0y0)) )
		)
	)
}

#' Calculates multiple confusion matrix stats from mariginal probabilities in dplyr friendly manner
#'
#' The purpose of this is to make it possible to calculate accuracy stats in a DPLYR table (including dbplyr). Typically
#' this will come from some sort of threshold based classification task where the classification output
#' is a probability and the prediction class is a binary outcome.
#' 
#' @param df a dataframe containing one observation per row & p_x1y1, p_x0y1, p_x1y0, and p_x0y0 columns (see probabilitiesFromConfusionMatrix)
#' @return the datatable with additional columns for true_pos_rate / true_neg_rate / etc...
#' @import dplyr
#' @export
calculateConfusionMatrixStats = function(df) {
	return(df %>% mutate(
		
		true_pos_rate = p_x1y1/p_x1,
		true_neg_rate = p_x0y0/p_x0,
		false_pos_rate = p_x0y1/p_x0,
		false_neg_rate = p_x1y0/p_x1,
		
		neg_pred_value = p_x0y0/p_y0,
		pos_pred_value = p_x1y1/p_y1,
		
		specificity = true_neg_rate,
		sensitivity = true_pos_rate,
				
		precision = pos_pred_value,
		recall = true_pos_rate,
		
		accuracy = p_x1y1+p_x0y0,
		f1 = 2*precision*recall/(precision+recall),
		mcc = (p_x1y1*p_x0y0 - p_x0y1*p_x1y0) / sqrt(p_x1*p_x0*p_y1*py0),
		informedness = true_pos_rate+true_neg_rate-1,
		youdens_j = informedness

	))
}








