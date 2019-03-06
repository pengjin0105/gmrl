#' Jacobian matrix calculation
#'
#'
#' @param data_order a ordered dataset by observed time
#' @param z_mat covariate matrix
#' @param beta0 parameter coefficients
#' @param model specify additive or proportional MRL model
#' @return jacobian matrix
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example


J_mat = function(data_order,z_mat,beta0,model='exp'){
  col = ncol(z_mat)
  h = matrix(NA,nrow = col, ncol = col)
  for(u in 1:col){
    for(be in 1:col){
      if(model=='exp'){h[u,be] = ub_exp(data_order,z_mat,beta0,u,be)}
      if(model=='add'){h[u,be] = ub_add(data_order,z_mat,beta0,u,be)}
    }
  }
  return(h)
}
