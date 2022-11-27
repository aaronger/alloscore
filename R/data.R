#' Data from Table 2 in ZXH for the normal distribution example
#'
#' A data frame with 17 rows and 9 columns:
#' \tabular{ll}{
#' 	\strong{Product} \tab	Product number \cr
#' 	\strong{v}		 \tab	cost of revenue loss per unit not stocked of product \emph{i}\cr
#' 	\strong{h}		 \tab	cost incurred per unit leftover of product \emph{i}	\cr
#' 	\strong{c}		 \tab	cost per unit of product \emph{i}\cr
#'	\strong{mu}		 \tab	mean of demand for product \emph{i}\cr
#'	\strong{sigma}	 \tab	standard deviation of demand for product \emph{i}\cr
#'	\strong{q_zxh}	 \tab	\emph{(v-c)/(h+v)} quantile of demand distribution as given by ZXH\cr
#'	\strong{GIM}	 \tab	solution of Abdel-Malek and Montanari (2005a)\cr
#'	\strong{Opt}	 \tab	solution of ZXH
#' }
#' @source B. Zhang, X. Xu, and Z. Hua, “A binary solution method for the multiproduct newsboy problem with budget constraint,” 
#' Int. J. Prod. Econ., vol. 117, no. 1, pp. 136–141, 2009.
"zxh_tab2"

#' Data from Table 3 in ZXH for the beta distribution example
#'
#' A data frame with 6 rows and 11 columns:
#' \tabular{ll}{
#' 	\strong{Product} \tab	Product number \cr
#' 	\strong{v}		 \tab	cost of revenue loss per unit not stocked of product \emph{i}\cr
#' 	\strong{h}		 \tab	cost incurred per unit leftover of product \emph{i}	\cr
#' 	\strong{c}		 \tab	cost per unit of product \emph{i}\cr
#'	\strong{x_min}	 \tab	min of demand distribution support\cr
#'	\strong{x_max}	 \tab	max of demand distribution support\cr
#'	\strong{Balpha}	 \tab	alpha (i.e. shape1) parameter of demand distribution\cr
#'	\strong{Bbeta}	 \tab	beta (i.e. shape2) parameter of demand distribution\cr
#'	\strong{GIM}	 \tab	solution of Abdel-Malek and Montanari (2005a)\cr
#'	\strong{Opt}	 \tab	solution of ZXH
#' }
#' @source B. Zhang, X. Xu, and Z. Hua, “A binary solution method for the multiproduct newsboy problem with budget constraint,” 
#' Int. J. Prod. Econ., vol. 117, no. 1, pp. 136–141, 2009.
"zxh_tab3"