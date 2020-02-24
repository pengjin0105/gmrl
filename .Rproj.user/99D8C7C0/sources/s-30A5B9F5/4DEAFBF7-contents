#' newton method
#'
#'
#' @param data_order a ordered dataset by observed time
#' @param z_mat covariate matrix
#' @param beta_ini given initial parameters
#' @param model specify additive or proportional MRL model
#' @param iter_max iteration maximum
#' @param tol tolerance
#' @return estimated value
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example




newton = function(data_order,z_mat,beta_ini,model='exp',iter_max=500,tol=0.0001){
  if(model=='exp'){
    for(j in 1:iter_max){
      newbeta = beta_ini - solve(J_mat(data_order,z_mat,beta_ini,model = 'exp'))%*%score_exp(data_order,z_mat,beta_ini)
      if(prod(abs(newbeta-beta_ini)<tol)==1){break}
      j=j+1
      beta_ini = newbeta
    }
    return(newbeta)
  }

  if(model=='add'){
    for(j in 1:iter_max){
      newbeta = beta_ini - solve(J_mat(data_order,z_mat,beta_ini,model = 'add'))%*%t(score_add(data_order,z_mat,beta_ini))
      if(prod(abs(newbeta-beta_ini)<tol)==1){break}
      j=j+1
      beta_ini = newbeta
    }
    return(newbeta)
  }
}
