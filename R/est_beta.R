#' main function
#'
#'
#' @param data_order a ordered dataset by observed time
#' @param z_mat covariate matrix
#' @param model specify model type
#' @return estimates and standard error of parameters
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example



est_beta = function(data_order,z_mat,model='exp'){
  library(survival)
  if(model=='exp'){
    tmp_fit = coxph(Surv(obs_time,delta)~z_mat,data=data_order)
    beta_ini = -summary(tmp_fit)$coef[,1]
    est = newton(data_order,z_mat,beta_ini,model='exp')
    se = vpart_exp(data_order,z_mat,est)
    return(list(estimate=est,
                se=se))
  }

  if(model=='add'){
    tmp_fit = lm(obs_time~z_mat,data=data_order)
    beta_ini = summary(tmp_fit)$coef[-1,1]
    est = newton(data_order,z_mat,beta_ini,model='add')
    se = vpart_add(data_order,z_mat,est)
    return(list(estimate=est,
                se=se))
  }
}
