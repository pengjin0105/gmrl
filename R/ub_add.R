#' Function to calculate derivative
#'
#'
#' @param data_order a ordered dataset by observed time
#' @param z_mat covariate matrix
#' @param beta0 parameter coefficients
#' @param u uth number of parameters
#' @param be bth number of parameters
#' @return derivative of dudb
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example


ub_add = function(data_order,z_mat,beta0,u,be){
  n = nrow(data_order)
  sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
  sn[is.na(sn)] = 0
  dm_tmp = -sn*data_order$pi*z_mat[,be]*data_order$delta/rev(cumsum(rev(data_order$pi)))
  dm_tmp[is.na(dm_tmp)] = 0
  dm = rev(cumsum(rev(dm_tmp)))/sn
  dm[is.na(dm)] = 0
  Z_bar = (apply(apply(apply(data_order$pi*z_mat,2,rev),2,cumsum),2,rev)/rev(cumsum(rev(data_order$pi))))[,u]
  Z_bar[is.na(Z_bar)] = 0
  ub_tmp = data_order$pi*(z_mat[,u]-Z_bar)*(dm+z_mat[,be])*data_order$delta
  return(sum(ub_tmp)/sqrt(n))
}
