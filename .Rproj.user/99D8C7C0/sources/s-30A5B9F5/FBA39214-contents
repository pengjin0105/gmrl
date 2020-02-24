#' standard error of estimated parameters
#'
#'
#' @param data_order a ordered dataset by observed time
#' @param z_mat covariate matrix
#' @param beta0 parameter coefficients
#' @return SEs of estimator
#' @author Peng Jin and Mengling Liu
#' @references
#' @export
#' @example



vpart_exp = function(data_order,z_mat,beta0){
  n = nrow(data_order)
  col = ncol(z_mat)
  bn = apply(apply(apply(data_order$pi*exp(-z_mat%*%beta0),2,rev),2,cumsum),2,rev)/
    rev(cumsum(rev(data_order$pi)))
  bn[is.na(bn)] = 0

  sn = exp(-cumsum(data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))))
  sn[is.na(sn)] = 0
  hat_tmp = sn[-n]*bn[-1]
  hat_m0 = rev(cumsum(rev(c(hat_tmp,0)*
                            diff(c(data_order$obs_time,data_order$obs_time[n])))))/sn
  hat_m0[is.na(hat_m0)] = 0

  z_bar = apply(apply(apply(data_order$pi*z_mat,2,rev),2,cumsum),2,rev)/rev(cumsum(rev(data_order$pi)))
  z_bar[is.na(z_bar)] = 0
  z_tild_tmp1 = data_order$delta*sn
  z_tild = apply(matrix(rep(z_tild_tmp1,col),ncol=col)*data_order$pi*(z_mat-z_bar),2,cumsum)*
    matrix(rep(sn/rev(cumsum(rev(data_order$pi))),col),ncol=col)
  z_tild[is.na(z_tild)] = 0

  A = matrix(0,col,col)
  v1 = matrix(0,col,col)
  v2 = matrix(0,col,col)
  v3 = matrix(0,col,col)
  vpart2_1 = matrix(0,col,col)
  vpart2_2 = matrix(0,1,col)
  for(k in 1:n){
    Z = z_mat[k,]

    A_tmp1 = matrix(rep(Z,each=n),n) - z_bar
    A_tmp2 = data_order$pi[k]*as.numeric(exp(-Z%*%beta0))*diff(c(0,data_order$obs_time))*
      (data_order$obs_time<=data_order$obs_time[k])
    A = A + t(A_tmp1)%*%(A_tmp1*A_tmp2)/n

    v1_tmp1 = matrix(rep(Z,each=n),n) - z_bar - z_tild
    v1_tmp2 = A_tmp2*hat_m0
    v1 = v1 + t(v1_tmp1)%*%(v1_tmp1*v1_tmp2)/n

    v2_tmp1 = data_order$pi[k]*data_order$delta/rev(cumsum(rev(data_order$pi)))*
      (data_order$obs_time<=data_order$obs_time[k])*hat_m0*hat_m0
    v2_tmp1[is.na(v2_tmp1)] = 0
    v2 = v2 + t(v1_tmp1)%*%(v1_tmp1*v2_tmp1)/n

    v3_tmp1 = data_order$pi[k]*(data_order$obs_time<=data_order$obs_time[k])*hat_m0*
      bn*diff(c(0,data_order$obs_time))
    v3 = v3 + t(v1_tmp1)%*%(v1_tmp1*matrix(rep(v3_tmp1,col),ncol=col))/n

    vpart2_tmp1 = (data_order$obs_time<=data_order$obs_time[k])*
      (as.numeric(exp(-Z%*%beta0))*diff(c(0,data_order$obs_time))+
         hat_m0*data_order$pi*data_order$delta/rev(cumsum(rev(data_order$pi)))-
         bn*diff(c(0,data_order$obs_time)))
    vpart2_tmp1[is.na(vpart2_tmp1)] = 0
    vpart2_tmp2 = apply(v1_tmp1*matrix(rep(vpart2_tmp1,col),ncol=col),2,sum)
    vpart2_1 = vpart2_1 + (data_order$pi[k]-data_order$delta[k])*vpart2_tmp2%*%t(vpart2_tmp2)/n
    vpart2_2 = vpart2_2 + (data_order$pi[k]-data_order$delta[k])*vpart2_tmp2/n
  }
  V1 = v1+v2-v3
  estp = sum(data_order$sub)/n
  V2 = (1-estp)/estp*(vpart2_1-t(vpart2_2)%*%vpart2_2)
  V = V1+V2
  return(sqrt(diag(solve(A)%*%V%*%solve(A)/n)))
}
