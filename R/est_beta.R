#' main function
#'
#'
#' @param time observed time variable in survival setting
#' @param delta censoring indicator (1=case and 0=censored subject)
#' @param vars covariate matrix
#' @param wt weight (unselected control has weight 0)
#' @param model specify model type (exp=proportional MRL model, add=additive MRL model)
#' @param se specify whether need calculating standard error or not (users can also use bootstrap or perturbation approach)
#' @return estimates and standard error of parameters
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example



est_beta = function(time,delta,vars,wt=1,model='exp',se.cal=FALSE){
  data = data.frame(obs_time=time,delta=delta,pi=wt,vars)
  data_order = data[order(data$obs_time),]
  z_mat = data_order[,4:ncol(data_order)]
  if(model=='exp'){
    ## estimate initial value from the Cox PH model
    tmp_fit = survival::coxph(Surv(obs_time,delta)~z_mat,data=data_order)
    beta_ini = -summary(tmp_fit)$coef[,1]
    est = newton(data_order,z_mat,beta_ini,model=model)
    if(se.cal==TRUE){
      se = vpart_exp(data_order,z_mat,est)
      return(list(estimate=est,
                  se=se))
    } else {
      return(list(estimate=est))
    }
  }

  if(model=='add'){
    ## estimate inital value from a linear regression
    tmp_fit = stats::lm(obs_time~z_mat,data=data_order)
    beta_ini = summary(tmp_fit)$coef[-1,1]
    est = newton(data_order,z_mat,beta_ini,model=model)
    se = vpart_add(data_order,z_mat,est)
    if(se.cal==TRUE){
      se = vpart_add(data_order,z_mat,est)
      return(list(estimate=est,
                  se=se))
    } else {
      return(list(estimate=est))
    }
  }
}
