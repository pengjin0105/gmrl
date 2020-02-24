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

vpart_add = function(data_order,z_mat,beta0){
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
  z_tild_tmp1 = data_order$delta*sn
  z_tild = apply(matrix(rep(z_tild_tmp1,col),ncol=col)*data_order$pi*(z_mat-z_bar),2,cumsum)*
    matrix(rep(sn/rev(cumsum(rev(data_order$pi))),2),ncol=2)
  z_tild[is.na(z_tild)] = 0

  # start estiamte A and Sigma
  A = matrix(0,col,col)
  v = matrix(0,col,col)
  vpart1 = matrix(0,col,col)
  vpart2 = matrix(0,1,col)
  vpart2_1 = matrix(0,col,col)
  vpart2_2 = matrix(0,1,col)
  for(k in 1:n){
    Z = z_mat[k,]
    A_tmp1 = matrix(rep(Z,each=n),n) - z_bar
    A = A + t(A_tmp1)%*%(A_tmp1*(data_order$obs_time==data_order$obs_time[k])*data_order$delta)*data_order$pi[k]/n
    sum_tmp = data_order$pi*(hat_m0+z_mat%*%beta0)*data_order$delta/rev(cumsum(rev(data_order$pi)))
    v_tmp1 = matrix(rep(Z,each=n),n) - z_bar - z_tild
    v_tmp2 = data_order$pi[k]*(hat_m0+as.numeric(z_mat[k,]%*%beta0))*
      sum_tmp*(data_order$obs_time<=data_order$obs_time[k])
    v_tmp2[is.na(v_tmp2)] = 0
    v = v + t(v_tmp1)%*%(v_tmp1*matrix(rep(v_tmp2,col),ncol=col))/n

    vpart1_tmp1 = (hat_m0+as.numeric(z_mat[k,]%*%beta0))*(data_order$obs_time==data_order$obs_time[k])*
      data_order$delta[k]
    vpart1_tmp2 = sum_tmp*(data_order$obs_time<=data_order$obs_time[k])
    vpart1_tmp3 = vpart1_tmp1-vpart1_tmp2
    vpart1_tmp3[is.na(vpart1_tmp3)] = 0
    vpart1_tmp4 = t(vpart1_tmp3)%*%v_tmp1
    vpart1 = vpart1 + (data_order$pi[k]-data_order$delta[k])*t(vpart1_tmp4)%*%(vpart1_tmp4)/n
    vpart2 = vpart2 + (data_order$pi[k]-data_order$delta[k])*vpart1_tmp4/n

  }
  estp = sum(data_order$sub)/n
  v2 = (1-estp)/estp*(vpart1-t(vpart2)%*%vpart2)
  V = v+v2
  return(sqrt(diag(solve(A)%*%V%*%solve(A)/n)))
}
