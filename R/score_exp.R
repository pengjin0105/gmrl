#' Score function for exponential MRL models
#'
#'
#' @param data_order a ordered dataset by observed time
#' @param z_mat covariate matrix
#' @param beta0 parameter coefficients
#' @return score for each parameter
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example


score_exp = function(data_order,z_mat,beta0){
  n = nrow(data_order)
  col = ncol(z_mat)
  sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
  sn[is.na(sn)] = 0
  bn = apply(apply(apply(data_order$pi*exp(-z_mat%*%beta0),2,rev),2,cumsum),2,rev)/
    rev(cumsum(rev(data_order$pi)))
  bn[is.na(bn)] = 0
  hat_tmp = sn[-n]*bn[-1]
  hat_m0 = rev(cumsum(rev(c(hat_tmp,0)*diff(c(data_order$obs_time,data_order$obs_time[n])))))/sn
  z_bar = apply(apply(apply(data_order$pi*z_mat,2,rev),2,cumsum),2,rev)/rev(cumsum(rev(data_order$pi)))
  exp_part = exp(-z_mat%*%beta0)
  zpart = diff(c(0,data_order$obs_time))*z_bar
  part1 = hat_m0*data_order$delta*(z_mat-z_bar)
  part1[is.na(part1)] = 0
  part2 = (z_mat*data_order$obs_time-apply(zpart,2,cumsum))
  part2 = part2*c(exp_part)
  part2[is.na(part2)] = 0
  return(colMeans(data_order$pi*(part1-part2)))
}
