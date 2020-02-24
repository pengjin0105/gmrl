#' Score function for additive MRL models
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


score_add = function(data_order,z_mat,beta0){
  n = nrow(data_order)
  col = ncol(z_mat)
  sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
  sn[is.na(sn)] = 0
  part1 = c(rev(cumsum(rev(sn[-1]*diff(data_order$obs_time)))),0)
  part2_tmp = sn*data_order$pi*z_mat%*%beta0*data_order$delta/rev(cumsum(rev(data_order$pi)))
  part2_tmp[is.na(part2_tmp)] = 0
  part2 = rev(cumsum(rev(part2_tmp)))
  hat_m0 = (part1-part2)/sn
  hat_m0[is.na(hat_m0)] = 0
  z_bar = apply(apply(apply(data_order$pi*z_mat,2,rev),2,cumsum),2,rev)/rev(cumsum(rev(data_order$pi)))
  z_bar[is.na(z_bar)] = 0
  u_tmp1 = data_order$pi*(hat_m0+z_mat%*%beta0)*data_order$delta
  u_tmp = t(u_tmp1)%*%(z_mat-z_bar)
  return(u_tmp/sqrt(n))
}
