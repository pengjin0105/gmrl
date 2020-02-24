#' Function to calculate derivative
#'
#'
#' @param data_order a ordered dataset by observed time
#' @param z_mat covariate matrix
#' @param beta0 parameter coefficients
#' @param u uth number of parameters
#' @param be bth number of parameters
#' @return element of Jacobian matrix
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example



ub_exp = function(data_order,z_mat,beta0,u,be){
  n = nrow(data_order)
  sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
  sn[is.na(sn)] = 0
  bn = apply(apply(apply(data_order$pi*exp(-z_mat%*%beta0)*(-z_mat[,be]),2,rev),2,cumsum),2,rev)/
    rev(cumsum(rev(data_order$pi)))
  bn[is.na(bn)] = 0
  dmb_tmp = sn[-n]*bn[-1]
  dmb = rev(cumsum(rev(c(dmb_tmp,0)*diff(c(data_order$obs_time,data_order$obs_time[n])))))/sn
  z_bar = (apply(apply(apply(data_order$pi*z_mat,2,rev),2,cumsum),2,rev)/rev(cumsum(rev(data_order$pi))))[,u]
  exp_part = exp(-z_mat%*%beta0)*(-z_mat[,be])
  zpart = diff(c(0,data_order$obs_time))*z_bar
  part1 = dmb*data_order$delta*(z_mat[,u]-z_bar)
  part1[is.na(part1)] = 0
  part2 = (z_mat[,u]*data_order$obs_time-cumsum(zpart))*exp_part
  part2[is.na(part2)] = 0
  return(mean(data_order$pi*(part1-part2)))
}
